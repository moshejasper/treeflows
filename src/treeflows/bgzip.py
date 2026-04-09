import io
import os
import struct
import zlib
from pathlib import Path


class BgzfWriter:
    """
    Minimal BGZF writer.

    - Writes valid BGZF blocks
    - Accepts bytes or str
    - Appends the standard BGZF EOF block on close
    - Does NOT build a .gzi index
    """

    # Conservative block payload size. BGZF allows up to 65536 uncompressed bytes
    # per block, but this slightly smaller value avoids edge-case oversize blocks.
    _MAX_UNCOMPRESSED_BLOCK = 65280

    # Standard 28-byte BGZF EOF marker from the HTS spec.
    _BGZF_EOF = bytes.fromhex(
        "1f8b08040000000000ff0600424302001b0003000000000000000000"
    )

    def __init__(self, filename, mode="wb", compresslevel=6, encoding="utf-8"):
        if "w" not in mode:
            raise ValueError("BgzfWriter only supports write modes.")
        if "b" not in mode:
            # allow "w" as shorthand, but file is always opened binary internally
            mode = mode.replace("t", "").replace("w", "wb")

        self.filename = filename
        self.compresslevel = compresslevel
        self.encoding = encoding
        self._buffer = bytearray()
        self._closed = False

        if isinstance(filename, (str, os.PathLike, Path)):
            self._fh = open(filename, "wb")
            self._owns_handle = True
        else:
            self._fh = filename
            self._owns_handle = False

    def write(self, data):
        if self._closed:
            raise ValueError("Cannot write to closed BgzfWriter.")

        if isinstance(data, str):
            data = data.encode(self.encoding)
        elif not isinstance(data, (bytes, bytearray, memoryview)):
            raise TypeError("data must be str or bytes-like")

        self._buffer.extend(data)

        while len(self._buffer) >= self._MAX_UNCOMPRESSED_BLOCK:
            block = bytes(self._buffer[:self._MAX_UNCOMPRESSED_BLOCK])
            del self._buffer[:self._MAX_UNCOMPRESSED_BLOCK]
            self._write_bgzf_block(block)

    def flush(self):
        if self._closed:
            return
        while self._buffer:
            block = bytes(self._buffer[:self._MAX_UNCOMPRESSED_BLOCK])
            del self._buffer[:self._MAX_UNCOMPRESSED_BLOCK]
            self._write_bgzf_block(block)
        self._fh.flush()

    def close(self):
        if self._closed:
            return
        self.flush()
        self._fh.write(self._BGZF_EOF)
        self._fh.flush()
        if self._owns_handle:
            self._fh.close()
        self._closed = True

    def __enter__(self):
        return me if (me := self) else self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def _write_bgzf_block(self, block):
        if len(block) > 65536:
            raise ValueError("A BGZF block cannot exceed 65536 uncompressed bytes.")

        # Raw DEFLATE stream for gzip payload
        compressor = zlib.compressobj(
            self.compresslevel,
            zlib.DEFLATED,
            -15,   # raw DEFLATE, not zlib wrapper
        )
        cdata = compressor.compress(block) + compressor.flush()

        # BGZF block = 18-byte header + compressed payload + 8-byte trailer
        total_block_size = 18 + len(cdata) + 8
        if total_block_size > 65536:
            raise RuntimeError(
                f"Compressed BGZF block too large ({total_block_size} bytes). "
                "Reduce _MAX_UNCOMPRESSED_BLOCK."
            )

        # gzip header with F.EXTRA set, XLEN=6
        header = bytearray()
        header.extend(b"\x1f\x8b")                  # ID1, ID2
        header.append(8)                           # CM = DEFLATE
        header.append(4)                           # FLG = F.EXTRA
        header.extend(struct.pack("<I", 0))        # MTIME
        header.append(0)                           # XFL
        header.append(255)                         # OS = unknown
        header.extend(struct.pack("<H", 6))        # XLEN = 6

        # BGZF extra subfield:
        # SI1='B', SI2='C', SLEN=2, BSIZE=total_block_size-1
        header.extend(b"BC")
        header.extend(struct.pack("<H", 2))
        header.extend(struct.pack("<H", total_block_size - 1))

        crc32 = zlib.crc32(block) & 0xFFFFFFFF
        isize = len(block) & 0xFFFFFFFF
        trailer = struct.pack("<II", crc32, isize)

        self._fh.write(header)
        self._fh.write(cdata)
        self._fh.write(trailer)
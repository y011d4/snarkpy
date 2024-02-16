from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass
class SectionPosSize:
    pos: int
    size: int


class Reader:
    def __init__(self, filename: Path) -> None:
        self._fp = filename.open("rb")

    def __del__(self) -> None:
        self._fp.close()

    def read_uint32(self) -> int:
        return self.read_uintn(4)

    def read_uint64(self) -> int:
        return self.read_uintn(8)

    def read_uintn(self, n: int) -> int:
        return int.from_bytes(self._fp.read(n), "little")

    def read(self, n: int) -> bytes:
        return self._fp.read(n)

    def get_pos(self) -> int:
        return self._fp.tell()

    def set_pos(self, n: int) -> int:
        return self._fp.seek(n)

    def skip(self, n: int) -> None:
        self._fp.seek(n, 1)

    def move_to_section(self, section: SectionPosSize) -> None:
        self.set_pos(section.pos)

    def check_section_end(self, section: SectionPosSize) -> bool:
        return self._fp.tell() == section.pos + section.size


class Sections:
    def __init__(self) -> None:
        self._sections: dict[int, list[SectionPosSize]] = defaultdict(list)

    def add_section(self, n: int, pos: int, size: int) -> None:
        self._sections[n].append(SectionPosSize(pos=pos, size=size))

    def get_section(self, n: int) -> SectionPosSize:
        assert len(self._sections[n]) == 1
        return self._sections[n][0]

    def get_sections(self, n: int) -> Sequence[SectionPosSize]:
        return self._sections[n]

    @classmethod
    def from_reader(cls, reader: Reader, file_type: bytes) -> "Sections":
        _type = reader.read(4)
        assert _type == file_type
        version = reader.read_uint32()
        assert version in [1, 2]
        n_sections = reader.read_uint32()
        sections = Sections()
        for _ in range(n_sections):
            ht = reader.read_uint32()
            hl = reader.read_uint64()
            sections.add_section(ht, reader.get_pos(), hl)
            reader.skip(hl)
        return sections

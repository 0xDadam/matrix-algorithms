from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, SupportsFloat


@dataclass(frozen=True, slots=True)
class CounterFloat:
    value: float
    additions: ClassVar[int] = 0
    subtractions: ClassVar[int] = 0
    multiplications: ClassVar[int] = 0
    divisions: ClassVar[int] = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", float(self.value))

    @classmethod
    def reset_counter(cls) -> None:
        cls.additions = 0
        cls.subtractions = 0
        cls.multiplications = 0
        cls.divisions = 0

    @classmethod
    def operation_counts(cls) -> dict[str, int]:
        return {
            "additions": cls.additions,
            "subtractions": cls.subtractions,
            "multiplications": cls.multiplications,
            "divisions": cls.divisions,
            "total": cls.additions
            + cls.subtractions
            + cls.multiplications
            + cls.divisions,
        }

    def _other(self, other: SupportsFloat) -> float:
        return float(other)

    def __float__(self) -> float:
        return self.value

    def __add__(self, other: SupportsFloat) -> CounterFloat:
        type(self).additions += 1
        return type(self)(self.value + self._other(other))

    def __radd__(self, other: SupportsFloat) -> CounterFloat:
        return self + other

    def __sub__(self, other: SupportsFloat) -> CounterFloat:
        type(self).subtractions += 1
        return type(self)(self.value - self._other(other))

    def __rsub__(self, other: SupportsFloat) -> CounterFloat:
        type(self).subtractions += 1
        return type(self)(self._other(other) - self.value)

    def __mul__(self, other: SupportsFloat) -> CounterFloat:
        type(self).multiplications += 1
        return type(self)(self.value * self._other(other))

    def __rmul__(self, other: SupportsFloat) -> CounterFloat:
        return self * other

    def __truediv__(self, other: SupportsFloat) -> CounterFloat:
        type(self).divisions += 1
        return type(self)(self.value / self._other(other))

    def __rtruediv__(self, other: SupportsFloat) -> CounterFloat:
        type(self).divisions += 1
        return type(self)(self._other(other) / self.value)

    def __repr__(self) -> str:
        return f"CounterFloat({self.value!r})"

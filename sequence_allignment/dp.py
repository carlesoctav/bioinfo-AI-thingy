from dataclasses import dataclass, field
from typing import Literal, Union


@dataclass
class PairAllignmentWithDPLinearGap:
    allignment_level: Literal["global", "local"]
    char_list: str = "AGCT"
    allign_score: int = 1
    mismatch_score: int = -1
    gap_score: int = -1
    contructed: bool = False
    table: Union[list[list[int]], str] = field(default_factory=list)
    solutions: list[tuple[str, str]] = field(default_factory=list)

    def construct_table(
        self,
        str_a: str,
        str_b: str,
    ):
        self.check_char_list(str_a, str_b)
        self.str_a = str_a
        self.str_b = str_b
        self.table = [[0 for _ in range(len(str_b) + 1)] for _ in range(len(str_a) + 1)]

        for i in range(1, len(str_a) + 1):
            self.table[i][0] = self.table[i - 1][0] + self.gap_score

        for j in range(1, len(str_b) + 1):
            self.table[0][j] = self.table[0][j - 1] + self.gap_score

        for i in range(1, len(str_a) + 1):
            for j in range(1, len(str_b) + 1):
                score = None
                if str_a[i - 1] == str_b[j - 1]:
                    score = self.allign_score
                else:
                    score = self.mismatch_score

                if self.allignment_level == "global":
                    self.table[i][j] = max(
                        self.table[i - 1][j - 1] + score,
                        self.table[i - 1][j] + self.gap_score,
                        self.table[i][j - 1] + self.gap_score,
                    )
                elif self.allignment_level == "local":
                    self.table[i][j] = max(
                        self.table[i - 1][j - 1] + score,
                        self.table[i - 1][j] + self.gap_score,
                        self.table[i][j - 1] + self.gap_score,
                        0,
                    )

        self.contructed = True

    def construct_allignment(self):
        if not self.contructed:
            raise ValueError("Table not constructed yet")

        start_i = len(self.str_a)
        start_j = len(self.str_b)
        queue = [(start_i, start_j, "", "", str(self.str_a), str(self.str_b))]

        while queue:
            print(f"==>> queue: {queue}")
            i, j, sol_a, sol_b, str_a, str_b = queue.pop(0)

            score = None

            if str_a == "" or str_b == "":
                score = -1e9
            elif str_a[-1] == str_b[-1]:
                score = self.allign_score
            else:
                score = self.mismatch_score

            if (
                i - 1 >= 0
                and j - 1 >= 0
                and self.table[i - 1][j - 1] + score == self.table[i][j]
            ):
                queue.append(
                    (
                        i - 1,
                        j - 1,
                        str_a[-1] + sol_a,
                        str_b[-1] + sol_b,
                        str_a[:-1],
                        str_b[:-1],
                    )
                )

            if (
                i - 1 >= 0
                and j >= 0
                and self.table[i - 1][j] + self.gap_score == self.table[i][j]
            ):
                queue.append(
                    (i - 1, j, str_a[-1] + sol_a, "-" + sol_b, str_a[:-1], str_b)
                )

            if (
                i >= 0
                and j - 1 >= 0
                and self.table[i][j - 1] + self.gap_score == self.table[i][j]
            ):
                queue.append(
                    (i, j - 1, "-" + sol_a, str_b[-1] + sol_b, str_a, str_b[:-1])
                )

            if str_a == "" and str_b == "":
                self.solutions.append((sol_a, sol_b))

    def check_char_list(
        self,
        str_a: str,
        str_b: str,
    ):
        if any([char not in self.char_list for char in str_a]):
            raise ValueError(
                f"Invalid character in str_a. Valid characters are {self.char_list}"
            )

        if any([char not in self.char_list for char in str_b]):
            raise ValueError(
                f"Invalid character in str_b. Valid characters are {self.char_list}"
            )

    def get_table(
        self,
    ) -> Union[list[list[int]], str]:
        if self.contructed:
            return self.table
        else:
            return "Table not constructed yet"

    def get_solutions(
        self,
    ) -> list[tuple[str, str]]:
        if self.contructed and self.solutions:
            return self.solutions
        else:
            return "Table or solutions not constructed yet"


if __name__ == "__main__":
    alligner = PairAllignmentWithDPLinearGap(
        allignment_level="global",
        char_list="AGCT",
        gap_score=-2,
        allign_score=1,
        mismatch_score=-1,
    )
    alligner.construct_table("CTTAGA", "GTAA")
    print(alligner.table)
    alligner.construct_allignment()
    print(alligner.solutions)

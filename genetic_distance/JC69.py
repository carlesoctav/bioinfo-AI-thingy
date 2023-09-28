import numpy as np


class JC69:
    def __init__(self, mu: float):
        self.mu = mu
        self.Q = np.array(
            [
                [-3 / 4 * self.mu, 1 / 4 * self.mu, 1 / 4 * self.mu, 1 / 4 * self.mu],
                [1 / 4 * self.mu, -3 / 4 * self.mu, 1 / 4 * self.mu, 1 / 4 * self.mu],
                [1 / 4 * self.mu, 1 / 4 * self.mu, -3 / 4 * self.mu, 1 / 4 * self.mu],
                [1 / 4 * self.mu, 1 / 4 * self.mu, 1 / 4 * self.mu, -3 / 4 * self.mu],
            ]
        )
        self.seq_list = []
        self.seq_name = []
        self.d_distance = None
        self.p_distance = None

    def get_prob_at_t(self, t: float):
        return np.linalg.matrix_power(np.exp(self.Q * t), t)

    def input(self, file: str):
        self.file_name = file.split(".")[0]
        with open(file, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if i == 0:
                    self.num_of_seq, self.seq_len = line.split()
                else:
                    name, seq = line.split()
                    self.seq_name.append(name)
                    self.seq_list.append(seq)

    def compute_p_distance(self):
        p_distance = np.zeros((len(self.seq_list), len(self.seq_list)))
        for i in range(len(self.seq_list)):
            for j in range(i + 1, len(self.seq_list)):
                p_distance[j][i] = self.compute_p_distance_between_two_seq(
                    self.seq_list[j], self.seq_list[i]
                )
        index = np.array(self.seq_name).reshape(-1, 1)

        p_distance_for_print = np.round(p_distance, 4)
        np.savetxt(
            f"{self.file_name}_p_distance.txt",
            np.column_stack((index, p_distance_for_print)),
            fmt="%s",
            delimiter="\t",
        )
        self.p_distance = p_distance
        return p_distance

    def compute_p_distance_between_two_seq(self, seq_a: str, seq_b: str):
        dissimilarity = 0
        for seq_1, seq_2 in zip(seq_a, seq_b):
            if seq_1 != seq_2:
                dissimilarity += 1
        return dissimilarity / len(seq_a)

    def compute_d_distance_exact(
        self,
    ):
        """
        since JC69 is a simple model, we can compute the d_distance directly using a given formula
        """
        p_distance = self.compute_p_distance()
        d_distance = np.zeros((len(self.seq_list), len(self.seq_list)))
        for i in range(len(self.seq_list)):
            for j in range(i + 1, len(self.seq_list)):
                d_distance[j][i] = -3 / 4 * np.log(1 - 4 / 3 * p_distance[j][i])
        index = np.array(self.seq_name).reshape(-1, 1)
        d_distance_for_print = np.round(d_distance, 4)
        np.savetxt(
            f"JC69_{self.file_name}_d_distance.txt",
            np.column_stack((index, d_distance_for_print)),
            fmt="%s",
            delimiter="\t",
        )
        self.d_distance = d_distance
        return d_distance

    def compute_difference(self):
        if self.d_distance is None:
            self.compute_d_distance_exact()
        if self.p_distance is None:
            self.compute_p_distance()

        difference = (self.d_distance - self.p_distance) / self.d_distance
        difference = np.nanmean(difference)
        print(f"==>> difference(in %): {difference* 100}")
        return difference


if __name__ == "__main__":
    model = JC69(0.1)
    model.input("primates.txt")
    model.compute_difference()

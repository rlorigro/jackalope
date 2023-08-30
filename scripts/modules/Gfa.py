

class GfaSequence:
    def __init__(self, tokens):
        self.name = tokens[1]
        self.sequence = tokens[2]

    def __len__(self):
        return len(self.sequence)


class GfaEdge:
    def __init__(self, tokens):
        self.name_a = tokens[1]
        self.reversal_a = self.get_reversal(tokens[2])

        self.name_b = tokens[3]
        self.reversal_b = self.get_reversal(tokens[4])

    def __str__(self):
        return self.name_a + ('-' if self.reversal_a else '+') + self.name_b + ('-' if self.reversal_b else '+')

    @staticmethod
    def get_reversal(token):
        if token == '-':
            return True
        elif token == '+':
            return False
        else:
            raise Exception("ERROR: unrecognized reversal token: " + token)


def iterate_gfa_nodes(gfa_path):
    with open(gfa_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            data = line.strip().split('\t')

            if data[0] == "S":
                yield GfaSequence(data)


def iterate_gfa_edges(gfa_path):
    with open(gfa_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            data = line.strip().split('\t')

            if data[0] == "L":
                yield GfaEdge(data)


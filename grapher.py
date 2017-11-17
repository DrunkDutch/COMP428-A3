import random
import sys
MAX_VALUE = 100
DISCO_CHANCE = 50
DISCO_RATE = 100


def gen_graph(inte):
    num_nodes = int(inte)
    with open("input.txt", 'w') as output:
        y = 0
        x = 0

        for y in range(num_nodes):
            for x in range(num_nodes):
                if x == y:
                    output.write("{}\t".format(0))
                    continue

                disconnect_chance = random.randint(0, 100) % DISCO_RATE
                if disconnect_chance > DISCO_CHANCE:
                    output.write("1\t")
                else:
                    output.write("inf\t")
            output.write("\n")


if __name__ == "__main__":
    gen_graph(sys.argv[1])

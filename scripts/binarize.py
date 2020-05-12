def to_binary(n, digits=8):
    l = []
    for _ in range(digits):
        l.append(n % 2)
        n = int(n / 2)
    return l


def make_dist_arr(n, inputs):
    outputs = to_binary(n)
    dist = []
    for i, inp in enumerate(inputs):
        dist.append(inp + str(outputs[i]))
    return dist

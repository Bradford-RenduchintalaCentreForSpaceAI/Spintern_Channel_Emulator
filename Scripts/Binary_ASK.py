import random

import matplotlib.pyplot as plt
import numpy as np


# We can include what's called "type hints" when defining function parameters.
# This way whoever uses the function will know to supply the correct kind of
# variable, and you wouldn't need to use things like `int(length)`
def generate_2ask(t, wt, num_samples: int, length: int, b1, step):
    """Please include some basic documentation with all your function

    Emphasis should be on what the parameters are, and restrictions or edge
    cases users need to be aware of. When describing the function, talk about
    *what* it does, rather than *how*. The *how* should be obvious from the
    code itself.

    Args:
        t (_type_): _description_
        wt (_type_): _description_
        num_samples (_type_): _description_
        length (_type_): _description_
        b1 (_type_): _description_
        step (_type_): _description_

    Returns:
        _type_: _description_
    """

    c = np.sin(t * wt)

    b1_mapped = []

    # This calculation only happens once, so no need to have it in the loops
    new_var = int((num_samples/length)/step)

    # Does this nested loop add multiple identical copies of `b1[q]` to
    # `b1_mapped`?
    for q in range(length):  # `range()` automatically starts from 0
        for i in range(new_var):
            b1_mapped.append(b1[q])

    s = []

    # If you use NumPy arrays there should be a better, quicker way to do
    # element-wise multiplication - `np.multiply` I think
    for i in range(len(t)):
        s.append(b1_mapped[i] * c[i])

    return s, b1_mapped


# Using this construct to turn this script into a reuseable module. This way
# if we just run the script we get the test case you have developed. However,
# if we import it into a different script or into a Jupyter notebook, the code
# below will not run and we can use the `generate_2ask` function on its own.
if __name__ == "__main__":

    num_samples = 10000

    # This is a bit of a guess - generally a good idea to use descriptive names
    length = num_samples / 8  # you can also use integer division, `//`
    step = 0.01

    b1 = random.choices([0, 1], weights=(50, 50), k=int(length))

    graph_scaling_factor = 5

    # In Python, `^` denotes XOR. Use `**` or `np.pow()` instead
    wt = (1 * 10**9) * np.pi

    t = np.arange(0, num_samples, step)

    [s, b1_mapped] = generate_2ask(t, wt, num_samples, int(length), b1, step)

    fig1, (sub1, sub2) = plt.subplots(1, 2, figsize=(10, 10))

    # Why does you plot discard the first element and start at index `1`?
    sub1.plot(
        [t[i] for i in range(1, int(num_samples / graph_scaling_factor))],
        [b1_mapped[i] for i in range(1, int(num_samples / graph_scaling_factor))]
    )

    sub1.plot(
        [t[i] for i in range(1, int(num_samples / graph_scaling_factor))],
        [s[i] for i in range(1, int(num_samples / graph_scaling_factor))],
        linestyle='dashed'
    )

    sub2.plot(np.fft.fft(b1_mapped))

    sub2.plot(np.fft.fft(s))

    sub2.legend(["Bit train", "Signal", "Carrier"])

    sub2.set(ylim=(0, 10000))

    plt.show()

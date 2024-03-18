import random
import string


def random_at(string):

    e = 0

    while e < len(string):
        if string[e : e + 5] == "[C@H]":
            first = string[:e]
            new_string = string[e:]
            string = new_string.replace(
                "[C@H]", random.sample(["[C@H*]", "[C@@H*]"], 1)[0], 1
            )
            string = first + string
            e = e + 5
        elif string[e : e + 4] == "[C@]":
            first = string[:e]
            new_string = string[e:]
            string = new_string.replace(
                "[C@]", random.sample(["[C@*]", "[C@@*]"], 1)[0], 1
            )
            string = first + string
            e = e + 4
        elif string[e : e + 6] == "[C@@H]":
            first = string[:e]
            new_string = string[e:]
            string = new_string.replace(
                "[C@@H]", random.sample(["[C@H*]", "[C@@H*]"], 1)[0], 1
            )
            string = first + string
            e = e + 6
        elif string[e : e + 5] == "[C@@]":
            first = string[:e]
            new_string = string[e:]
            string = new_string.replace(
                "[C@@]", random.sample(["[C@*]", "[C@@*]"], 1)[0], 1
            )
            string = first + string
            e = e + 4
        else:
            e = e + 1

    string = string.replace("*", "")
    return string

import rearrangement
import time

now = time.gmtime()
displaytime = time.strftime("%d_%m_%Y__%H_%M", now)
print(displaytime)

f = open("test2", "r")

g = open("smiles_good_order", "w")

ligne = f.readline()

while ligne != "":
    a = rearrangement.rearrangement(ligne)
    # print a
    if a != None:
        g.write(a[0])

    ligne = f.readline()

now2 = time.gmtime()
displaytime2 = time.strftime("%d_%m_%Y__%H_%M", now2)
print(displaytime2)
f.close()
g.close()

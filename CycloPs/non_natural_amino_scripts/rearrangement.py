import string
import random
import re
import create_liste
import random_at
from rdkit import Chem


def rearrangement(string):
    """
    Rearranges SMILES strings to N-term-[middle]-C-term
    Input string should be of the form SMILES\tName, such as below.
    [*]N(C)CC(=O)[O-]	ZINC00057605
    """
    m = open("f.log", "a")
    zinc = []
    zinc.extend(re.findall("\t(ZINC\d+)\n", string))
    smiles = []
    smiles.extend(re.findall("\A(.+)\t", string))
    smiles = smiles[0]

    # Keep N at the start and the carboxylate at the end
    a = smiles.replace("C(=O)[O-]", "")

    b = re.search("N\(.+\)\[\*\]", a)

    if "[*]N" in a:
        c = a.replace("[*]N", "")
    elif "N[*]" in a:
        c = a.replace("N[*]", "")
    elif b != None:
        chaine = a[::-1]
        e = 0
        while e < len(chaine):
            if chaine[e : e + 4] == "]*[)":
                i = e + 4
                s = ""
                while chaine[i : i + 2] != "(N":
                    s = s + chaine[i]
                    i = i + 1
                s = ")" + s + "("
                s = s[::-1]

                c = re.sub("N\(.+\)\[\*\]", s, a)
                e = e + 1
            else:
                e = e + 1

    side_chain = c

    test = False
    counter = 0

    while test == False and counter < 10000:
        print((1))

        side_chain = c
        # times=random.randint(0,side_chain.count('('))

        # side_chain=side_chain.replace('(','',times)
        # side_chain=side_chain.replace(')','',times)

        liste = create_liste.create_liste(side_chain)

        random.shuffle(liste)
        # print liste

        random_list = liste

        final_string = "".join(random_list)

        final_string = random_at.random_at(final_string)

        final_string = "[*]N" + final_string + "C(=O)[O-]"

        check = [smiles, final_string]

        # print check, counter

        m.write(final_string + "\n")
        print((2))

        try:
            cans = []
            for smi in check:
                mol = Chem.MolFromSmiles(smi)
                print(mol)
                cans.append(Chem.MolToSmiles(mol, True))
            print((cans, 3))
            if cans[0] == cans[1]:
                print((4))
                # del cans
                test = True
            else:
                print((5))
                # del cans
                # del check
                test = False
        except:
            counter = counter + 1
            print((6))
            test = False

    if test == False:
        g = open("Others", "a")
        g.write(string)
        g.close()
        m.close()
    else:
        final_string = final_string + "\t" + zinc[0] + "\n"
        m.close()
        return final_string, counter

    # rearrangement('CN(c1ccc(cc1)C(=O)[O-])[*]\tZINC00000000\n')

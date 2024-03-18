def create_liste(string):
    liste = []
    dechet = ""
    compt = 0

    i = 0
    while i < len(string) - 1:
        if string[i] in "hHLBbcCnNoOfFMApPsSkKTvVZGiIRYyXwWDEUu()123456789":
            # case of atoms with one letter followed by a special character or a number
            if string[i + 1] in "[]()=#@":
                liste = liste + [string[i]]
                i = i + 1
            elif string[i + 1] in "1234567890":
                liste = liste + [string[i : i + 2]]
                i = i + 2
            # case of atoms with two letters
            elif string[i : i + 2] in ["Cl", "Br"]:
                liste = liste + [string[i : i + 2]]
                i = i + 2
            else:
                liste = liste + [string[i]]
                i = i + 1
        # case of special charater except [ or a number
        elif string[i] in "\/":
            i = i + 1
        # case in braquets : ions, chiral atoms, ...
        elif string[i] == "[":
            s = ""
            h = i
            lag = 0
            while string[h] != "]":
                s = s + string[h]
                lag = lag + 1
                h = h + 1
            s = s + "]"
            liste = liste + [s]
            i = i + lag

        # case of double or triple bounds
        elif string[i] in "#=":
            liste = liste + [string[i : i + 2]]
            i = i + 2

        # other case
        else:
            dechet = dechet + string[i]
            compt = compt + 1
            i = i + 1

    if string[len(string) - 1] in "hHLBbcCnNoOfFMApPsSkKTvVZGiIRYyXwWDEUu()123456789":
        liste = liste + [string[len(string) - 1]]

    return liste

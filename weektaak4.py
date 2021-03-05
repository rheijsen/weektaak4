from textwrap import wrap


def lees_inhoud(bestand_naam):
    headers_int = []
    seq_int = []
    headers_opp = []
    seq_opp = []
    seq = ""
    with open(bestand_naam) as fasta:
        for line in fasta:
            line = line.strip()
            if ">" in line:
                if seq != "":
                    seq_int.append(seq)
                    seq = ""
                headers_int.append(line)
            else:
                seq += line.strip()
                seq_int.append(seq)
    c = 0
    for i in headers_int:
        if "envelope" in i or "protein=env" in i:
            headers_opp.append(i)
            seq_opp.append(seq_int[c])
            del seq_int[c]
            del headers_int[c]
            c += 1
    return seq_int, seq_opp


def list_to_string(seq_int, seq_opp):
    string_int_seq = "".join(seq_int)
    string_opp_seq = "".join(seq_opp)
    return string_int_seq, string_opp_seq


def split_seq(string_int, string_opp):
    int_protein = wrap(string_int, 1)
    opp_protein = wrap(string_opp, 1)
    return int_protein, opp_protein


def count_amino_acids(int_protein, opp_protein):
    int_set = [[i, int_protein.count(i)] for i in set(int_protein)]
    opp_set = [[i, opp_protein.count(i)] for i in set(opp_protein)]
    return int_set, opp_set


def gegevens_int(int_set, length_int, string_int):
    aantal_c = string_int.count("C")
    print("Percentage Cysteine: ", round(100 * aantal_c/ length_int, 2))
    aantal_w = string_int.count("W")
    print("Percentage Tryptofaan: ", round(100 * aantal_w / length_int, 2))
    aantal_hydrofoob = string_int.count("G") + string_int.count("A") + string_int.count("V") + string_int.count("L") + string_int.count("L") + string_int.count("I") + string_int.count("P") + string_int.count("F") + string_int.count("M") + string_int.count("W")
    print("Percentage hydrofoob: ", round(100 * aantal_hydrofoob / length_int, 2))
    aantal_hydrofiel = string_int.count("K") + string_int.count("R") + string_int.count("E") + string_int.count("D") + string_int.count("Q") + string_int.count("N")
    print("Percentage hydrofiel: ", round(100 * aantal_hydrofiel / length_int, 2))
    sorted_list = (sorted(int_set, key = lambda x: x[1]))
    print("Meest voorkomende aminozuren:")
    print(sorted_list[19][0], round(100 * sorted_list[19][1] / length_int, 2))
    print(sorted_list[18][0], round(100 * sorted_list[18][1] / length_int, 2))
    print(sorted_list[17][0], round(100 * sorted_list[17][1] / length_int, 2))
    print("Minst voorkomende aminozuren:")
    print(sorted_list[0][0], round(100 * sorted_list[0][1] / length_int, 2))
    print(sorted_list[1][0], round(100 * sorted_list[1][1] / length_int, 2))
    print(sorted_list[2][0], round(100 * sorted_list[2][1] / length_int, 2))
    print("-" * 80)


def gegevens_opp(opp_set, length_opp, string_opp):
    print("Oppevlakte proteinen: \n")
    aantal_c = string_opp.count("C")
    print("Percentage Cysteine: ", round(100 * aantal_c/ length_opp, 2))
    aantal_w = string_opp.count("W")
    print("Percentage Tryptofaan: ", round(100 * aantal_w / length_opp, 2))
    aantal_hydrofoob = string_opp.count("G") + string_opp.count("A") + string_opp.count("V") + string_opp.count("L") + string_opp.count("L") + string_opp.count("I") + string_opp.count("P") + string_opp.count("F") + string_opp.count("M") + string_opp.count("W")
    print("Percentage hydrofoob: ", round(100 * aantal_hydrofoob / length_opp, 2))
    aantal_hydrofiel = string_opp.count("K") + string_opp.count("R") + string_opp.count("E") + string_opp.count("D") + string_opp.count("Q") + string_opp.count("N")
    print("Percentage hydrofiel: ", round(100 * aantal_hydrofiel / length_opp, 2))
    sorted_list = (sorted(opp_set, key = lambda x: x[1]))
    print("Meest voorkomende aminozuren:")
    print(sorted_list[18][0], round(100 * sorted_list[18][1] / length_opp, 2))
    print(sorted_list[17][0], round(100 * sorted_list[17][1] / length_opp, 2))
    print(sorted_list[16][0], round(100 * sorted_list[16][1] / length_opp, 2))
    print("Minst voorkomende aminozuren:")
    print(sorted_list[0][0], round(100 * sorted_list[0][1] / length_opp, 2))
    print(sorted_list[1][0], round(100 * sorted_list[1][1] / length_opp, 2))
    print(sorted_list[2][0], round(100 * sorted_list[2][1] / length_opp, 2))


if __name__ == '__main__':
    bestand = "HIV-2 Protein.txt"
    seq_int, seq_opp = lees_inhoud(bestand)
    string_int_seq, string_opp_seq = list_to_string(seq_int, seq_opp)
    int_protein, opp_protein = split_seq(string_int_seq, string_opp_seq)
    int_set, opp_set = count_amino_acids(int_protein, opp_protein)
    gegevens_int(int_set, len(string_int_seq), string_int_seq)
    gegevens_opp(opp_set, len(string_opp_seq), string_opp_seq)

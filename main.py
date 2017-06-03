import random
import math
import cmath
import sys

epsilon = math.pow(10, -9)
h = math.pow(10, -6)
k_max = 10000


class Pereche:
    a = None
    b = None

    def __init__(self, a, b):
        self.a = a
        self.b = b


def main():
    pb1()
    pb2()


def pb1():
    p = list()
    parametri_polinom = open('polinom.txt', "rt").read().split(";")

    coeficienti = list()
    for t in range(len(parametri_polinom) - 1, -1, -1):
        coeficienti.append(float(parametri_polinom[t]))

    R = (abs(coeficienti[len(coeficienti) - 1]) + max(coeficienti)) / abs(coeficienti[len(coeficienti) - 1])

    print("Interval radacini: [" + str(-R) + "," + str(R) + "]")
    print("Metoda lui Laguerre: ")

    p = list()
    for i in range(0, len(coeficienti)):
        # print(coeficienti[i], sep=', ')
        p.append(complex(coeficienti[i]))
    # print(p)

    radacini = gaseste_radacini(p, R)
    # print(roots)
    size_radacini = len(radacini)
    radacini_unice = list()

    for i in range(0, size_radacini):
        if (abs(radacini[i].real) < epsilon):
            radacini[i] = radacini[i] - complex(radacini[i].real, 0)
        if (abs(radacini[i].imag) < epsilon):
            radacini[i] = radacini[i] - complex(0, radacini[i].imag)
        radacini.append(complex(radacini[i].real, -radacini[i].imag))

    radacini.sort(key=lambda x: x.real)
    radacini_unice.append(radacini[0])

    for i in range(1, len(radacini)):
        already_exists = False
        for ii in range(0, len(radacini_unice)):
            if (abs(radacini_unice[ii].real - radacini[i].real) <= epsilon) and (
                        abs(radacini_unice[ii].imag - radacini[i].imag) <= epsilon):
                already_exists = True
                break
        if already_exists == False:
            radacini_unice.append(radacini[i])

    print("")

    fisier_radacini = open('radacini.txt', "wt")

    print("Radacinile sunt: ")
    for i in range(0, len(radacini_unice)):
        print(str(radacini_unice[i].real) + ", " + str(radacini_unice[i].imag) + "i")
        fisier_radacini.write(str(radacini_unice[i].real) + ", " + str(radacini_unice[i].imag) + "i" + "\n")


def pb2():
    parametri_polinom = open('polinom.txt', "rt").read().split(";")
    p = list()
    for t in range(0, len(parametri_polinom)):
        p.append(float(parametri_polinom[t]))

    R = (abs(p[0]) + max(p)) / abs(p[0])

    solutii1 = list()
    solutii_unice1 = list()

    print("")
    print("Formula de aproximare G1")

    for count_incercari1 in range(0, 5):
        x_prec = random.randint(-int(R), int(R))
        x = x_prec + 1
        k = 0

        while True:
            fx_derivat = calculeaza_derivata_G1(p, x)
            fx_prec_derivat = calculeaza_derivata_G1(p, x_prec)
            numitor_delta_x = (fx_derivat - fx_prec_derivat)
            if numitor_delta_x >= -epsilon and numitor_delta_x <= epsilon:
                delta_x = math.pow(10, -5)
            else:
                delta_x = ((x - x_prec) * fx_derivat) / numitor_delta_x
            x_prec = x
            x = x - delta_x
            k = k + 1
            if abs(delta_x) < epsilon or k > k_max or abs(delta_x) > math.pow(10, 8):
                break
        if abs(delta_x) < epsilon:
            solutii1.append(x)
            if len(solutii_unice1) == 0:
                solutii_unice1.append(x)
            else:
                already_exists = False
                for contor in range(0, len(solutii_unice1)):
                    if (abs(solutii_unice1[contor] - x) <= epsilon):
                        already_exists = True
                        break
                if already_exists == False:
                    solutii_unice1.append(x)
        else:
            print("divergenta")
            count_incercari1 = count_incercari1 - 1

    for s in range(0, len(solutii_unice1)):
        derivata_secunda = calculeaza_derivata_secunda(p, solutii_unice1[s])
        if derivata_secunda > 0:
            print("Punctul critic " + str(solutii_unice1[s]) + " este punct de minim")
        else:
            print("Punctul critic " + str(solutii_unice1[s]) + " NU este punct de minim")

    solutii2 = list()
    solutii_unice2 = list()

    print("")
    print("Formula de aproximare G2")

    count_incercari2 = 0
    counter = 0
    while count_incercari2 < 5:
        x_prec = random.randint(-int(R), int(R))
        while True:
            x = random.randint(-int(R), int(R))
            if abs(x - x_prec) >= epsilon:
                break
        k = 0

        while True:
            fx_derivat = calculeaza_derivata_g2(p, x)
            fx_prec_derivat = calculeaza_derivata_g2(p, x_prec)
            numitor_delta_x = (fx_derivat - fx_prec_derivat)
            if numitor_delta_x >= -epsilon and numitor_delta_x <= epsilon:
                delta_x = math.pow(10, -5)
            else:
                delta_x = ((x - x_prec) * fx_derivat) / numitor_delta_x

            x_prec = x
            x = x - delta_x
            k = k + 1
            if abs(delta_x) < epsilon or k > k_max or abs(delta_x) > math.pow(10, 8):
                break
        if abs(delta_x) < epsilon:
            solutii2.append(x)
            if len(solutii_unice2) == 0:
                solutii_unice2.append(x)
            else:
                already_exists = False
                for contor in range(0, len(solutii_unice2)):
                    if (abs(solutii_unice2[contor] - x) <= epsilon):
                        already_exists = True
                        break
                if already_exists == False:
                    solutii_unice2.append(x)
        else:
            count_incercari2 = count_incercari2 - 1

        count_incercari2 = count_incercari2 + 1
        counter = counter + 1
        if counter >= 50:
            print("Divergenta")
            break

    for s in range(0, len(solutii_unice2)):
        derivata_secunda = calculeaza_derivata_secunda(p, solutii_unice2[s])
        if derivata_secunda > 0:
            print("Punctul critic " + str(solutii_unice2[s]) + " este punct de minim")
        else:
            print("Punctul critic " + str(solutii_unice2[s]) + " NU este punct de minim")


def gaseste_radacini(p, R):
    q = p
    w = q[:]

    rez = list()
    for i in range(0, 10):
        q = w[:]
        while len(q) > 2:
            z = complex(random.randint(-int(R), int(R)), random.randint(-int(R), int(R)))
            z = gaseste_radacina(q, z)
            z = gaseste_radacina(p, z)
            q = horner(q, z).a
            q.pop()
            rez.append(z)
        rez.append(-q[0] / q[1])
    return rez


def gaseste_radacina(p0, x):
    n = len(p0) - 1
    p1 = deriveaza(p0)
    p2 = deriveaza(p1)
    for step in range(0, 10000):
        y0 = eval(p0, x)
        if compara(y0, 0) == 0:
            break
        G = eval(p1, x) / y0
        H = G * G - eval(p2, x) - y0
        R = cmath.sqrt(complex(n - 1) * (H * complex(n) - G * G))
        D1 = G + R
        D2 = G - R
        if (compara(D1, D2)) > 0:
            toDivide = D1
        else:
            toDivide = D2
        a = complex(n) / toDivide
        x = x - a
        if (compara(a, 0) == 0):
            break

    return x


def deriveaza(p):
    n = len(p)
    rez = list()
    for i in range(1, n):
        rez.append(p[i] * complex(i))
    return rez


def eval(p, x):
    return horner(p, x).b;


def horner(a, x0):
    n = len(a)
    b = [0] * n

    for i in range(n - 1, 0, -1):
        to_add = 0
        if i < (n - 1):
            to_add = b[i] * x0
        b[i - 1] = a[i] + to_add

    return Pereche(b, a[0] + b[0] * x0)


def compara(x, y):
    diff = abs(x) - abs(y)
    if diff < -epsilon:
        return -1
    else:
        if diff > epsilon:
            return 1
        else:
            return 0


def get_lista_derivare(a):
    lista_derivare = [0] * (len(a) - 1)
    for i in range(0, len(lista_derivare)):
        lista_derivare[i] = a[i] * (a.length - 1 - i);
    return lista_derivare


def calculeaza_horner_polinomial(a, x):
    b = a[0];
    for i in range(1, len(a)):
        b = a[i] + b * x;
    return b


def calculeaza_derivata_G1(a, x):
    F1 = calculeaza_horner_polinomial(a, x)
    F2 = calculeaza_horner_polinomial(a, x - h)
    F3 = calculeaza_horner_polinomial(a, x - 2 * h)
    return (3 * F1 - 4 * F2 + F3) / h


def calculeaza_derivata_g2(a, x):
    F1 = calculeaza_horner_polinomial(a, x + 2 * h)
    F2 = calculeaza_horner_polinomial(a, x + h)
    F3 = calculeaza_horner_polinomial(a, x - h)
    F4 = calculeaza_horner_polinomial(a, x - 2 * h)
    return (-F1 + 8 * F2 - 8 * F3 + F4) / 12 * h


def calculeaza_derivata_secunda(a, x):
    F1 = calculeaza_horner_polinomial(a, x + 2 * h)
    F2 = calculeaza_horner_polinomial(a, x + h)
    F3 = calculeaza_horner_polinomial(a, x)
    F4 = calculeaza_horner_polinomial(a, x - h)
    F5 = calculeaza_horner_polinomial(a, x - 2 * h)
    return (-F1 + 16 * F2 - 30 * F3 + 16 * F4 - F5) / 12 * h * h


main()

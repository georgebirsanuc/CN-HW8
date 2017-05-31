import random
import math
import cmath
import sys

epsilon = math.pow(10, -9)
h = math.pow(10, -6)
kMax = 10000


class Pair:
    first = None
    second = None

    def __init__(self, pFirst, pSecond):
        self.first = pFirst
        self.second = pSecond


def main():
    pb1()
    pb2()


def pb1():
    p = list()
    parametriPolinom = open('polinom.txt', "rt").read().split(";")

    coefficients = list()
    for t in range(len(parametriPolinom) - 1, -1, -1):
        coefficients.append(float(parametriPolinom[t]))

    R = (abs(coefficients[len(coefficients) - 1]) + max(coefficients)) / abs(coefficients[len(coefficients) - 1])

    print("Intervalul in care se gasesc radacinile polinomului este  [" + str(-R) + "," + str(R) + "]")
    print("Se aplica metoda lui Laguerre ...")

    p = list()
    for i in range(0, len(coefficients)):
        p.append(complex(coefficients[i]))
    print(p)

    roots = find_roots(p, R)
    print(roots)
    size_roots = len(roots)
    unique_roots = list()

    for i in range(0, size_roots):
        if (abs(roots[i].real) < epsilon):
            roots[i] = roots[i] - complex(roots[i].real, 0)
        if (abs(roots[i].imag) < epsilon):
            roots[i] = roots[i] - complex(0, roots[i].imag)
        roots.append(complex(roots[i].real, -roots[i].imag))

    roots.sort(key=lambda x: x.real)
    unique_roots.append(roots[0])

    for i in range(1, len(roots)):
        alreadyExists = False
        for contor in range(0, len(unique_roots)):
            if (abs(unique_roots[contor].real - roots[i].real) <= epsilon) and (
                        abs(unique_roots[contor].imag - roots[i].imag) <= epsilon):
                alreadyExists = True
                break
        if alreadyExists == False:
            unique_roots.append(roots[i])

    print("")

    polinoameFile = open('radacini.txt', "wt")

    print("Radacinile sunt: ")
    for i in range(0, len(unique_roots)):
        print(str(unique_roots[i].real) + " , " + str(unique_roots[i].imag) + "i")
        polinoameFile.write(str(unique_roots[i].real) + " , " + str(unique_roots[i].imag) + "i" + "\n")


def pb2():
    parametriPolinom = open('polinom.txt', "rt").read().split(";")
    p = list()
    for t in range(0, len(parametriPolinom)):
        p.append(float(parametriPolinom[t]))

    R = (abs(p[0]) + max(p)) / abs(p[0])

    solutii1 = list()
    solutiiUnice1 = list()

    print("")
    print("Formula de aproximare G1")

    for countIncercari1 in range(0, 5):
        xPrec = random.randint(-int(R), int(R))
        x = xPrec + 1
        k = 0

        while True:
            FdeXderivat = calculateDerivataG1(p, x)
            FdeXprecDerivat = calculateDerivataG1(p, xPrec)
            numitorDeltaX = (FdeXderivat - FdeXprecDerivat)
            if numitorDeltaX >= -epsilon and numitorDeltaX <= epsilon:
                deltaX = math.pow(10, -5)
            else:
                deltaX = ((x - xPrec) * FdeXderivat) / numitorDeltaX
            xPrec = x
            x = x - deltaX
            k = k + 1
            if abs(deltaX) < epsilon or k > kMax or abs(deltaX) > math.pow(10, 8):
                break
        if abs(deltaX) < epsilon:
            solutii1.append(x)
            if len(solutiiUnice1) == 0:
                solutiiUnice1.append(x)
            else:
                alreadyExists = False
                for contor in range(0, len(solutiiUnice1)):
                    if (abs(solutiiUnice1[contor] - x) <= epsilon):
                        alreadyExists = True
                        break
                if alreadyExists == False:
                    solutiiUnice1.append(x)
        else:
            print("divergenta")
            countIncercari1 = countIncercari1 - 1

    for r in range(0, len(solutiiUnice1)):
        derivataSecunda = calculateDerivateSecunda(p, solutiiUnice1[r])
        if derivataSecunda > 0:
            print("Punctul critic " + str(solutiiUnice1[r]) + " este punct de minim")
        else:
            print("Punctul critic " + str(solutiiUnice1[r]) + " NU este punct de minim")

    solutii2 = list()
    solutiiUnice2 = list()

    print("")
    print("Formula de aproximare G2")

    countIncercari2 = 0
    counter = 0
    while countIncercari2 < 5:

        xPrec = random.randint(-int(R), int(R))
        while True:
            x = random.randint(-int(R), int(R))
            if abs(x - xPrec) >= epsilon:
                break
        k = 0

        while True:
            FdeXderivat = calculateDerivataG2(p, x)
            FdeXprecDerivat = calculateDerivataG2(p, xPrec)
            numitorDeltaX = (FdeXderivat - FdeXprecDerivat)
            if numitorDeltaX >= -epsilon and numitorDeltaX <= epsilon:
                deltaX = math.pow(10, -5)
            else:
                deltaX = ((x - xPrec) * FdeXderivat) / numitorDeltaX

            xPrec = x
            x = x - deltaX
            k = k + 1
            if abs(deltaX) < epsilon or k > kMax or abs(deltaX) > math.pow(10, 8):
                break
        if abs(deltaX) < epsilon:
            solutii2.append(x)
            if len(solutiiUnice2) == 0:
                solutiiUnice2.append(x)
            else:
                alreadyExists = False
                for contor in range(0, len(solutiiUnice2)):
                    if (abs(solutiiUnice2[contor] - x) <= epsilon):
                        alreadyExists = True
                        break
                if alreadyExists == False:
                    solutiiUnice2.append(x)
        else:
            countIncercari2 = countIncercari2 - 1

        countIncercari2 = countIncercari2 + 1
        counter = counter + 1
        if counter >= 50:
            print("Divergenta")
            break

    for r in range(0, len(solutiiUnice2)):
        derivataSecunda = calculateDerivateSecunda(p, solutiiUnice2[r])
        if derivataSecunda > 0:
            print("Punctul critic " + str(solutiiUnice2[r]) + " este punct de minim")
        else:
            print("Punctul critic " + str(solutiiUnice2[r]) + " NU este punct de minim")


def find_roots(p, R):
    q = p
    w = q[:]

    res = list()
    for i in range(0, 10):
        q = w[:]
        while len(q) > 2:
            z = complex(random.randint(-int(R), int(R)), random.randint(-int(R), int(R)))
            z = find_a_root(q, z)
            z = find_a_root(p, z)
            q = horner(q, z).first
            q.pop()
            res.append(z)
        res.append(-q[0] / q[1])
    return res


def find_a_root(p0, x):
    n = len(p0) - 1
    p1 = derivare(p0)
    p2 = derivare(p1)
    for step in range(0, 10000):
        y0 = eval(p0, x)
        if compare(y0, 0) == 0:
            break
        G = eval(p1, x) / y0
        H = G * G - eval(p2, x) - y0
        R = cmath.sqrt(complex(n - 1) * (H * complex(n) - G * G))
        D1 = G + R
        D2 = G - R
        if (compare(D1, D2)) > 0:
            toDivide = D1
        else:
            toDivide = D2
        a = complex(n) / toDivide
        x = x - a
        if (compare(a, 0) == 0):
            break

    return x


def derivare(p):
    n = len(p)
    r = list()
    for i in range(1, n):
        r.append(p[i] * complex(i))
    return r


def eval(p, x):
    return horner(p, x).second;


def horner(a, x0):
    n = len(a)
    b = [0] * n

    for i in range(n - 1, 0, -1):
        toAdd = 0
        if i < (n - 1):
            toAdd = b[i] * x0
        b[i - 1] = a[i] + toAdd

    return Pair(b, a[0] + b[0] * x0)


def compare(x, y):
    diff = abs(x) - abs(y)
    if diff < -epsilon:
        return -1
    else:
        if diff > epsilon:
            return 1
        else:
            return 0


def getListDerivare(a):
    derivative = [0] * (len(a) - 1)
    for i in range(0, len(derivative)):
        derivative[i] = a[i] * (a.length - 1 - i);
    return derivative


def calculatePolynomialHorner(a, x):
    b = a[0];
    for i in range(1, len(a)):
        b = a[i] + b * x;
    return b


def calculateDerivataG1(a, x):
    F1 = calculatePolynomialHorner(a, x)
    F2 = calculatePolynomialHorner(a, x - h)
    F3 = calculatePolynomialHorner(a, x - 2 * h)
    return (3 * F1 - 4 * F2 + F3) / h


def calculateDerivataG2(a, x):
    F1 = calculatePolynomialHorner(a, x + 2 * h)
    F2 = calculatePolynomialHorner(a, x + h)
    F3 = calculatePolynomialHorner(a, x - h)
    F4 = calculatePolynomialHorner(a, x - 2 * h)
    return (-F1 + 8 * F2 - 8 * F3 + F4) / 12 * h


def calculateDerivateSecunda(a, x):
    F1 = calculatePolynomialHorner(a, x + 2 * h)
    F2 = calculatePolynomialHorner(a, x + h)
    F3 = calculatePolynomialHorner(a, x)
    F4 = calculatePolynomialHorner(a, x - h)
    F5 = calculatePolynomialHorner(a, x - 2 * h)
    return (-F1 + 16 * F2 - 30 * F3 + 16 * F4 - F5) / 12 * h * h


main()

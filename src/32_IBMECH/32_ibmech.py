'''
# change the default MNT224 to SS512 for the paper settings
# group = PairingGroup('MNT224', secparam=1024)
# group = PairingGroup('SS512', secparam=1024)
#
Shorter IBE and Signatures via Asymmetric Pairings
'''
from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.matrixops import *
from charm.toolbox.IBEnc import IBEnc
import copy
import time

class IBE_Chen12_z(IBEnc):

    def __init__(self, groupObj):
        IBEnc.__init__(self)
        global group, n
        group = groupObj
        n = 8

    def setup(self):
        start = time.time()
        g1 = group.random(G1)
        g2 = group.random(G2)
        alpha = group.random(ZR)
        eta = group.random(ZR)

        d1n, d2n, d3n, d4n, d5n, d6n,d7n, d8n = dnn_func()
        one = group.random(ZR)

        D1n, D2n, D3n, D4n = Dn_func(d1n, d2n,d3n, d4n, d5n, d6n, d7n,d8n, one)

        PP2 = (pair(g1, g2)) ** (alpha * one)
        PP3 = (pair(g1, g2)) ** (eta * one)

        gd1n = gd_func(g1, d1n)
        gd2n = gd_func(g1, d2n)
        gd3n = gd_func(g1, d3n)
        gd4n = gd_func(g1, d4n)

        pk = {'PP2': PP2, 'PP3': PP3, 'gd1n': gd1n, 'gd2n': gd2n}

        msk = {'alpha': alpha, 'eta': eta, 'g1': g1, 'g2': g2,
               'gd3n': gd3n, 'gd4n': gd4n, 'D1n': D1n, 'D2n': D2n, 'D3n': D3n, 'D4n': D4n,
               }

        end = time.time()
        rt = end - start
        return pk, msk, rt

    def skgen(self, msk, ID):
        start = time.time()

        _ID = group.hash(ID)
        r = group.random(ZR)
        ek_idn = ek_id_func(msk, _ID, r)
        k = {'ek_idn': ek_idn}

        end = time.time()
        rt = end - start
        return k, rt

    def rkgen(self, pk, msk, IDstar):
        start = time.time()

        _IDstar = group.hash(IDstar)
        s, s1, s2 = group.random(ZR, 3)
        sk_id1n = sk_id1n_func(msk, _IDstar, s, s1)
        sk_id2n = sk_id2n_func(msk, _IDstar, s, s2)
        sk_id3 =  pk['PP3'] ** s
        sk = {'sk_id1n': sk_id1n, 'sk_id2n': sk_id2n, 'sk_id3': sk_id3}

        end = time.time()
        rt = end - start
        return sk, rt

    def encrypt(self, pk, ek, IDstar, M):
        start = time.time()

        z = group.random(ZR)
        _IDstar = group.hash(IDstar)
        C0 = (pk['PP2'] ** z) * M
        C1n = C1n_func(pk, ek, z, _IDstar)
        ct = {'C0': C0, 'C1n': C1n}

        end = time.time()
        rt = end - start
        return ct, rt

    def decrypt(self, sk,  ID,  ct):
        start = time.time()
        _ID = group.hash(ID)
        a1 =a1_func(sk, _ID, ct)
        Mprime = ct['C0'] * sk['sk_id3'] / a1

        end = time.time()
        rt = end - start
        return Mprime, rt

def dn_func():
    res = []
    for i in range(n):
        res.append(group.random(ZR))
    return res

def dnn_func():
    d1n = dn_func()
    d2n = dn_func()
    d3n = dn_func()
    d4n = dn_func()
    d5n = dn_func()
    d6n = dn_func()
    d7n = dn_func()
    d8n = dn_func()

    return d1n, d2n, d3n, d4n, d5n, d6n, d7n, d8n


def Dn_func(d1n, d2n,d3n, d4n, d5n, d6n, d7n,d8n, one):
    d1n_, d2n_, d3n_, d4n_, d5n_, d6n_, d7n_, d8n_ = copy.copy(d1n), copy.copy(d2n), copy.copy(d3n), copy.copy(d4n), copy.copy(d5n), copy.copy(d6n), copy.copy(d7n), copy.copy(d8n)
    d1n_.append(group.init(ZR, 0))
    d2n_.append(group.init(ZR, 0))
    d3n_.append(group.init(ZR, 0))
    d4n_.append(group.init(ZR, 0))
    d5n_.append(group.init(ZR, 0))
    d6n_.append(group.init(ZR, 0))
    d7n_.append(group.init(ZR, 0))
    d8n_.append(group.init(ZR, 0))

    d1none, d2none, d3none, d4none = copy.copy(d1n), copy.copy(d2n), copy.copy(d3n), copy.copy(d4n)
    d1none.append(one)
    d2none.append(one)
    d3none.append(one)
    d4none.append(one)

    D1n = GaussEliminationinGroups([d1none, d2n_, d3n_, d4n_, d5n_, d6n_, d7n_, d8n_])
    D2n = GaussEliminationinGroups([d1n_, d2none, d3n_, d4n_, d5n_, d6n_, d7n_, d8n_])
    D3n = GaussEliminationinGroups([d1n_, d2n_, d3none, d4n_, d5n_, d6n_, d7n_, d8n_])
    D4n = GaussEliminationinGroups([d1n_, d2n_, d3n_, d4none, d5n_, d6n_, d7n_, d8n_])

    return D1n, D2n, D3n, D4n

def gd_func(g1, d1n):
    res = []
    for i in range(n):
        res.append(g1 ** d1n[i])
    return res


def ek_id_func(msk, _ID, r):
    res =[]
    for i in range(n):
        res.append(msk['gd3n'][i] ** (msk['eta'] + r * _ID) / (msk['gd4n'][i] ** r))
    return res


def sk_id1n_func(msk, _IDstar, s, s1):
    res = []
    for i in range(n):
        res.append(msk['g2'] ** (msk['D1n'][i] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D2n'][i] + s * msk['D3n'][i]))
    return res


def sk_id2n_func(msk, _IDstar, s, s2):
    res = []
    for i in range(n):
        res.append(msk['g2'] ** (s2 * (_IDstar * msk['D1n'][i] - msk['D2n'][i]) + s * msk['D4n'][i]))
    return res


def C1n_func(pk, ek, z, _IDstar):
    res = []
    for i in range(n):
        res.append((pk['gd1n'][i] ** z) * (pk['gd2n'][i] ** (z * _IDstar)) * (ek['ek_idn'][i]))
    return res


def a1_func(sk, _ID, ct):
    res = group.init(ZR, 1)
    for i in range(n):
        res = res * (pair(ct['C1n'][i], sk['sk_id1n'][i] * (sk['sk_id2n'][i] ** _ID)))
    return res


def main():
    group = PairingGroup('SS512', secparam=1024)
    ibe = IBE_Chen12_z(group)
    (pk, msk, setuptime) = ibe.setup()
    ID = 'user@email.com'
    IDstar = 'name@gmail.com'

    ek, skgentime = ibe.skgen(msk, ID)
    sk, rkgentime = ibe.rkgen(pk, msk, IDstar)
    msg = group.random(GT)
    print("ori message: ", msg)
    ct, enctime = ibe.encrypt(pk, ek, IDstar, msg)
    dec_msg, dectime = ibe.decrypt(sk, ID,  ct)
    print("dec message: ", dec_msg)
    print(dec_msg == msg)

    output_txt = './32_IBMECH.txt'
    with open(output_txt, 'w+', encoding='utf-8') as f:
        f.write("Seq SetupAveTime       ekgenAveTime       rkgenAveTime       encAveTime         decAveTime        " + '\n')
        f.write(str(n).zfill(2) + '  ' + str(setuptime) + ' ' + str(skgentime) + ' ' + str(rkgentime) + ' ' + str(enctime) + ' ' + str(dectime))
        f.write('\n')


if __name__ == '__main__':
    main()

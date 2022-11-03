# -*- coding: utf-8 -*-
"""
# message
# https://github.com/cygnusv/matchmaking-encryption/blob/master/ibme.py
# H0: {0, 1}* -> G0   H_prime()
# H1: {0, 1}* -> G1   H_prime()
# 'hello1234' -> '123456789'
# '123456789' -> 'hello1234'
# refer to VFPPBA

# from charm.toolbox.hash_module import Hash
# H = Hash(group)
# H.hashToZr()
# H.hashtoZn()
# H2: {0 ,1}* -> Zp  H.hashToZr()
# H3: GT -> Zp       H.hashtoZn()
# H4: phi = H.hashToZr(d1, d2, C1, C0, C2, C3)


"""
import random
from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEncMultiAuth import ABEncMultiAuth
import time
import numpy as np
from charm.toolbox.hash_module import Hash
from charm.core.math.integer import int2Bytes
import sys
import string
import random

sys.path.append('../')
from common.image import *
from common.msp import *
from common.id import *


class MJ18(ABEncMultiAuth):
    def __init__(self, groupObj, verbose=False):
        ABEncMultiAuth.__init__(self)
        global group, ahnipe, util, H, mask, id_len, l, str_len, l1

        group = groupObj
        util = MSP(group, verbose=False)
        H = Hash(group)
        mask = 'ed27dbfb02752e0e16bc4502d6c732bc5f1cc92ba19b2d93a4e95c597ca42753e93550b52f82b6c13fb8cc0c2fc64487'
        id_len = 2
        l = 50
        str_len = 48
        l1 = 15

    def setup_psbme(self, n):
        start = time.time()
        g = group.random(G1)
        h, u, v, w = group.random(G2), group.random(G2), group.random(G2), group.random(G2)
        alpha, beta, rou = group.random(ZR), group.random(ZR), group.random(ZR)
        g1 = g ** rou
        h0 = h ** rou
        h1 = h ** beta

        pp = {'g': g, 'g1': g1, 'u': u, 'v': v, 'w': w, 'h': h, 'h0': h0, 'h1': h1}
        msk = {"rou": rou, 'alpha': alpha}

        end = time.time()
        rt = end - start
        return pp, msk, rt

    def ekgen_psbme(self, msk, idstar):
        start = time.time()

        ekid = H1(idstar) ** msk['alpha']
        ek = {'ekid': ekid}

        end = time.time()
        rt = end - start
        return ek, rt

    def dkgen_psbme(self, msk, id):
        start = time.time()
        dk1 = H0(id) ** msk["rou"]
        dk2 = H0(id) ** msk["alpha"]
        dk3 = H0(id)
        dk = {'dk1': dk1, 'dk2': dk2, 'dk3': dk3}

        end = time.time()
        rt = end - start

        return dk, rt

    def enc_psbme(self, pp, ek, id_num):
        start = time.time()

        m = "123456789"
        s, d1, d2, sigma, tau = group.random(ZR), group.random(ZR), group.random(ZR), group.random(ZR), group.random(ZR)
        d1, d1cut = cutd_func(d1, 2)
        d2, d2cut = cutd_func(d2, 2)

        C0 = pp['h'] ** s
        C1 = pp['g'] ** s
        C2 = pp['h1'] ** tau

        C3temp = H3(d1, d2, C1, C0, C2)
        C3 = C3_func(C3temp, m)

        H = Hash(group)
        phi = H.hashToZr(C1, C0, C2, C3)
        C4 = ((pp['u'] ** phi) * (pp['v'] ** sigma) * (pp['w'])) ** s

        idj = idj_func(id_num)  # idj means S

        Uidtemp = Uid_func(pp, idj, s)
        Uid = cut_func(Uidtemp, 2)

        Vidtemp = Vid_func(ek, idj, C2)
        Vid = cut_func(Vidtemp, 2)

        an = fxgy_func(Uid, d1cut)
        bn = fxgy_func(Vid, d2cut)

        ct = {'sigma': sigma, 'C0': C0, 'C1': C1, 'C2': C2, 'C3': C3, 'C4': C4, 'an': an, 'bn': bn, 'idj': idj}
        end = time.time()
        rt = end - start
        return ct, m, rt

    def dec_psbme(self, pp, dk, idstar, ct):
        start = time.time()

        phi = H.hashToZr(ct['C1'], ct["C0"], ct["C2"], ct['C3'])

        a1 = pair(ct['C1'], (pp['u'] ** phi) * (pp['v'] ** ct['sigma']) * pp["w"])
        a2 = pair(pp['g'], ct['C4'])

        if a1 == a2:
            # print("pair correct verification")
            Uidtemp = H2(pair(ct['C0'], dk['dk1']))
            Uid = int(str(Uidtemp)[:2])
            Vidtemp = H2(pair(dk['dk3'], ct['C2']) * pair(dk['dk2'], H1(idstar)))
            Vid = int(str(Vidtemp)[:2])

            d1dec = ddec_func(ct['an'], Uid)
            d2dec = ddec_func(ct['bn'], Vid)

            a3 = str(ct["C3"][:str_len - l1])
            c3l = str(ct["C3"][str_len - l1:str_len])
            a4 = str(H3(d1dec, d2dec, ct['C1'], ct["C0"], ct['C2']))[:str_len - l1]
            c4l = str(H3(d1dec, d2dec, ct['C1'], ct["C0"], ct['C2']))[str_len - l1: str_len]
            if a3 == a4:
                # print("test dec successfully")
                dec_msg = int(c4l) ^ int(c3l)

        end = time.time()
        rt = end - start
        return dec_msg, rt


def H0(X):
    X = bytes([a ^ b for (a, b) in zip(X.encode(), bytes.fromhex(mask))])
    return group.hash(X, G1)


def H1(X):
    X = bytes([a ^ b for (a, b) in zip(X.encode(), bytes.fromhex(mask))])
    return group.hash(X, G2)


def H2(X):
    H = Hash(group)
    return int(H.hashToZr(X))


def H3(X1, X2, X3, X4, X5):
    a1, a2, a3, a4, a5 = int(str(X1)), int(str(X2)), int(H2(X3)), int(H2(X4)), int(H2(X5))
    a6 = str(a1 ^ a2 ^ a3 ^ a4 ^ a5)
    H = Hash(group)
    a7 = H0(a6)
    return H.hashToZn(a7)


def id_generate_func(len):
    res = ''.join(random.choices(string.ascii_uppercase +
                                 string.digits, k=len))
    return res


def idj_func(id_num):
    res = []
    for i in range(id_num):
        res.append(id_generate_func(id_len))
    return res


def fxgy_func(Vid, d1):
    fx = coeff_func(Vid)
    temp = fx[len(fx) - 1] + int(d1)
    fx[len(fx) - 1] = temp
    return fx


def cut_func(H2idjtemp, k):
    res = []
    for i in range(len(H2idjtemp)):
        res.append(int(str(H2idjtemp[i])[:k]))
    return res


def cutd_func(ele, k):
    d1cut = int(str(ele)[:k])
    d1 = group.init(ZR, int(str(d1cut)))
    return d1, d1cut


def coeff_func(H2idj):
    res = np.array([1, -H2idj[0]])
    for i in range(len(H2idj) - 1):
        a2 = np.array([1, -H2idj[i + 1]])
        res = np.convolve(res, a2)
    return res


def Uid_func(pp, idj, s):
    res = []
    for i in range(len(idj)):
        a1 = H0(idj[i]) ** s
        res.append(int(H2(pair(pp['h0'], a1))))
    return res


def Vid_func(ek, idj, C2):
    res = []
    for i in range(len(idj)):
        a1 = H0(idj[i])
        a2 = ek['ekid'] * C2
        res.append(int(H2(pair(a1, a2))))
    return res


def ddec_func(an, Uid):
    dd = 0
    for i in range(len(an)):
        dd = dd + an[i] * (Uid ** (len(an) - i - 1))
    res = int2elment(dd)
    return res


def int2elment(val):
    return group.init(ZR, int(val))


def C3_func(C3temp, m):
    n0 = str_len - l1
    C31 = str(C3temp)[0:n0]
    C32 = str(C3temp)[n0:str_len]
    C3 = C31 + str(int(C32) ^ int(m))
    return C3


def main():
    groupObj = PairingGroup('SS512')
    n_array = np.arange(2, 30, 1)
    # n_array = [10]
    output_txt = './31_PSBME.txt'
    ahnipe = MJ18(groupObj)

    with open(output_txt, 'w+', encoding='utf-8') as f:
        f.write(
            "Seq SetupAveTime       ekgenAveTime       EncAveTime         dkgenAVeTime       decAveTime        " + '\n')

        for i in range(len(n_array)):
            seq = 1
            sttot, ekgentot, enctot, dkgentot, dectot = 0.0, 0.0, 0.0, 0.0, 0.0
            for j in range(seq):
                n = n_array[i]
                idstar = id_generator(n)

                pp, msk, setuptime = ahnipe.setup_psbme(n)
                ek, ekgentime = ahnipe.ekgen_psbme(msk, idstar)
                ct, m, enctime = ahnipe.enc_psbme(pp, ek, n)
                dk, dkgentime = ahnipe.dkgen_psbme(msk, ct['idj'][0])
                rec_msg, dectime = ahnipe.dec_psbme(pp, dk, idstar, ct)

                print('\nn, seq:   ', n, j)
                print("m:       ", m)
                print("rec_msg: ", rec_msg)

                # m_inputkey = group.serialize(m).decode("utf-8")
                # m_inputkey = group.serialize(m).decode("utf-8")
                # m_outputkey = group.serialize(rec_msg).decode("utf-8")
                # encrypt(m_inputkey)
                # decrypt(m_outputkey)
                # image.encrypt(m_inputkey)
                # image.decrypt(m_outputkey)

                sttot, ekgentot, enctot, dkgentot, dectot = sttot + setuptime, ekgentot + ekgentime, enctot + enctime, dkgentot + ekgentime, dectot + dectime

            out0 = str(n).zfill(2)
            out1 = str(format(sttot / float(seq), '.16f'))
            out2 = str(format(ekgentot / float(seq), '.16f'))
            out3 = str(format(enctot / float(seq), '.16f'))
            out4 = str(format(dkgentot / float(seq), '.16f'))
            out5 = str(format(dectot / float(seq), '.16f'))
            f.write(out0 + '  ' + out1 + ' ' + out2 + ' ' + out3 + ' ' + out4 + ' ' + out5)
            f.write('\n')


if __name__ == "__main__":
    main()

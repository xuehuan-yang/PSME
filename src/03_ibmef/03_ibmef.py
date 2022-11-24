# -*- coding: utf-8 -*-
"""
# H0: {0, 1}* -> G0   H_prime()
# H1: {0, 1}* -> G1   H_prime()


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
from common.parse import *


class MJ18(ABEncMultiAuth):
    def __init__(self, groupObj, verbose=False):
        ABEncMultiAuth.__init__(self)
        global group, ahnipe, util, H, mask
        group = groupObj
        util = MSP(group, verbose=False)
        H = Hash(group)
        mask = 'ed27dbfb02752e0e16bc4502d6c732bc5f1cc92ba19b2d93a4e95c597ca42753e93550b52f82b6c13fb8cc0c2fc64487'

    def setup_psbme(self, n):
        start = time.time()
        g = group.random(G1)
        h = group.random(G2)
        alpha, y = group.random(ZR, 2)
        galpha = g ** alpha
        h1 = h ** y

        pp = {'g': g, 'galpha': galpha, 'h': h, 'h1': h1}
        msk = {'alpha': alpha, 'y': y}

        end = time.time()
        rt = end - start
        return pp, msk, rt

    def skgen_psbme(self, idstar):
        start = time.time()

        ek = {"idstar": idstar}

        end = time.time()
        rt = end - start
        return ek, rt

    def rkgen_psbme(self, pp, msk, id):
        start = time.time()

        H = Hash(group)
        rho = H.hashToZr(id)
        r = group.random(ZR)
        hrho = pp["h"] ** ((msk['y'] - r) / (msk['alpha'] - rho))
        dk = {'r': r, 'hrho': hrho}

        end = time.time()
        rt = end - start
        return dk, rt

    def enc_psbme(self, pp, ek, id):
        start = time.time()
        m = group.random(GT)
        s = group.random(ZR)
        x, enc_idstar = ext_func(pp, ek['idstar'])
        H = Hash(group)
        c1 = (pp['galpha'] * (pp['g'] ** (- H.hashToZr(id)))) ** s
        c2 = pair(pp['g'], pp['h']) ** s
        c3 = x
        c4 = m * (pair(pp['g'], pp['h1']) ** (-s)) * enc_idstar
        ct = {'c1': c1, "c2": c2, 'c3': c3, 'c4': c4}

        end = time.time()
        rt = end - start
        return ct, m, rt

    def dec_psbme(self, pp, dk, idstar, ct):
        start = time.time()

        a1 = pair(pp['g'], pp['h']) ** (H.hashToZr(idstar) * ct['c3'])
        dec_msg = (ct['c4'] * pair(ct['c1'], dk['hrho'])) * (ct['c2'] ** dk['r']) / a1

        end = time.time()
        rt = end - start
        return dec_msg, rt


def ext_func(pp, idstar):
    x = group.random(ZR)
    enc_idstar = pair(pp['g'], pp['h']) ** (x * H.hashToZr(idstar))
    return x, enc_idstar


def main():
    groupObj = PairingGroup('SS512')
    n_array = np.arange(5, 55, 5)
    n_array = np.insert(n_array, 0, 1)
    output_txt= output_func('33_ibmef')
    ahnipe = MJ18(groupObj)

    with open(output_txt, 'w+', encoding='utf-8') as f:
        f.write(
            "Seq SetupAveTime       ekgenAveTime       rkgenAveTime       encAVeTime         decAveTime         " + '\n')

        for i in range(len(n_array)):
            seq = seq_func()
            sttot, ekgentot, rkgentot, enctot, dectot = 0.0, 0.0, 0.0, 0.0, 0.0
            for j in range(seq):
                n = n_array[i]
                idstar = id_generator(n)
                id = id_generator(n)

                pp, msk, setuptime = ahnipe.setup_psbme(n)
                ek, ekgentime = ahnipe.skgen_psbme(idstar)
                dk, rkgentime = ahnipe.rkgen_psbme(pp, msk, id)
                ct, m, enctime = ahnipe.enc_psbme(pp, ek, id)
                rec_msg, dectime = ahnipe.dec_psbme(pp, dk, idstar, ct)

                print('\nn, seq:   ', n, j)
                print("m:        ", m)
                print("rec_msg1: ", rec_msg)

                m_inputkey = group.serialize(m).decode("utf-8")
                m_outputkey = group.serialize(rec_msg).decode("utf-8")
                encrypt(m_inputkey)
                decrypt(m_outputkey)

                sttot, ekgentot, rkgentot, enctot, dectot = sttot + setuptime, ekgentot + ekgentime, rkgentot + rkgentime, enctot + enctime, dectot + dectime

            out0 = str(n).zfill(2)
            out1 = str(format(sttot / float(seq), '.16f'))
            out2 = str(format(ekgentot / float(seq), '.16f'))
            out3 = str(format(rkgentot / float(seq), '.16f'))
            out4 = str(format(enctot / float(seq), '.16f'))
            out5 = str(format(dectot / float(seq), '.16f'))
            f.write(out0 + '  ' + out1 + ' ' + out2 + ' ' + out3 + ' ' + out4 + ' ' + out5)
            f.write('\n')


if __name__ == "__main__":
    main()

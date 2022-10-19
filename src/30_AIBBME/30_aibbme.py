# -*- coding: utf-8 -*-
"""
# https://github.com/LoneWalker/Implementation-CCGRID-CPABE/blob/621d1180ea650650e23d93e01b386ddb1a09537c/charm/uhibe_l12.py


# https://github.com/cygnusv/matchmaking-encryption/blob/master/ibme.py
# H0: {0, 1}* -> G1   H_prime()
# H1: {0, 1}* -> G0   H_prime()

# from charm.toolbox.hash_module import Hash
# H = Hash(group)
# H.hashToZr()
# H.hashtoZn()
# H2: {0 ,1}* -> Zp  H.hashToZr()
# H3: GT -> Zp       H.hashtoZn()

"""
import random
from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEncMultiAuth import ABEncMultiAuth
import time
import numpy as np
from charm.toolbox.hash_module import Hash
import sys
import string
import random

sys.path.append('../')
from common.image import *
from common.msp import *


class MJ18(ABEncMultiAuth):
    def __init__(self, groupObj, verbose=False):
        ABEncMultiAuth.__init__(self)
        global group, ahnipe, util, H, mask, id_len, id_num, l

        group = groupObj
        util = MSP(group, verbose=False)
        H = Hash(group)
        mask = 'ed27dbfb02752e0e16bc4502d6c732bc5f1cc92ba19b2d93a4e95c597ca42753e93550b52f82b6c13fb8cc0c2fc64487'
        id_len = 2
        id_num = 10
        l = 20

    def setup_aibbme(self):
        start = time.time()
        g, v = group.random(G1), group.random(G1)
        h = group.random(G2)
        r1n, r2n = random_array(l + 1), random_array(l + 1)
        t1, t2, beta1, beta2, alpha, rho = group.random(ZR, 6)
        b, tau = group.random(ZR, 2)

        rn = r_func(l + 1, b, r1n, r2n)
        t = t1 + b * t2
        beta = beta1 + b * beta2
        Rn = R_func(l + 1, g, rn)
        T = g ** t
        X = pair(g, h) ** beta

        vrho = v ** rho
        gb = g ** b
        hr1n = hr1_func(l + 1, h, r1n)
        hr2n = hr2_func(l + 1, h, r2n)
        ht1 = h ** t1
        ht2 = h ** t2
        gtaubeta = g ** (tau * beta)
        htaubeta1 = h ** (tau * beta1)
        htaubeta2 = h ** (tau * beta2)
        h1dtau = h ** (1 / tau)

        hbeta1 = h ** beta1
        hbeta2 = h ** beta2

        pp = {'v': v, 'vrho': vrho, 'g': g, 'gb': gb, 'R': Rn, 'T': T, 'X': X, 'h': h, 'hr1n': hr1n, 'hr2n': hr2n,
              'ht1': ht1, 'ht2': ht2, 'gtaubeta': gtaubeta, 'htaubeta1': htaubeta1, 'htaubeta2': htaubeta2,
              'h1dtau': h1dtau, 'r1ntest': r1n, 'r2ntest': r2n}
        msk = {'hbeta1': hbeta1, 'hbeta2': hbeta2, 'alpha': alpha, 'rho': rho}

        end = time.time()
        rt = end - start
        return pp, msk, rt

    def ekgen_aibbme(self, msk, idstar):
        start = time.time()

        ekid = H1(idstar) ** msk['alpha']
        ek = {'ekid': ekid, 'idstar': idstar}

        end = time.time()
        rt = end - start
        return ek, rt

    def dkgen_aibbme(self, pp, msk, id, ct):
        start = time.time()

        z = group.random(ZR)
        rtagsn = rtagsn_func(l)
        dk1 = H0(id) ** msk["rho"]
        dk2 = H0(id) ** msk["alpha"]
        dk3 = H0(id)
        dk4 = msk["hbeta1"] * (pp['ht1'] ** z)
        dk5 = msk["hbeta2"] * (pp['ht2'] ** z)
        dk6 = pp['h'] ** z
        dk7n = dk7_func(pp, rtagsn, id, z, ct)

        dk8n = dk8_func(pp, rtagsn, id, z, ct)
        dk = {'dk1': dk1, 'dk2': dk2, 'dk3': dk3, 'dk4': dk4, 'dk5': dk5, 'dk6': dk6, 'dk7n': dk7n, 'dk8n': dk8n,
              'rtagsn': rtagsn, 'z': z}

        end = time.time()
        rt = end - start

        return dk, rt

    def enc_aibbme(self, pp, ek):
        start = time.time()

        idj = idj_func(id_num)
        H2idjtemp = H2idj_func(idj)
        H2idj = cut_func(H2idjtemp, 2)
        ytemp = coeff_func(H2idj)
        ytest = yelement_func(ytemp)
        y = ytest
        print("enc y: ", y)

        # generate message
        m = group.random(GT)
        print("ori message: ", m)

        s, d2ori, ctag = group.random(ZR, 3)
        d2 = group.init(ZR, int(str(d2ori)[:5]))
        print("enc d2: ", d2)

        C0 = m * pp['X'] ** s
        C1 = pp["g"] ** s
        C2 = pp["gb"] ** s
        C3 = C3_func(pp, d2, s, ctag, ytest)
        C4 = pp['v'] ** s
        Vidtemp = Vid_func(pp, ek, idj, s)
        Vid = cut_func(Vidtemp, 2)
        print("enc Vid: ", Vid)

        C3prime = C3prime_func(pp, s, ctag, y)
        vidtemp = coeff_func(Vid)
        vidtemp1 = vidtemp
        bk = bk_func(vidtemp1, d2)

        ct = {'C0': C0, 'C1': C1, 'C2': C2, 'C3': C3, 'C4': C4, 'ctag': ctag, 'bk': bk, 'idj': idj, 'd2': d2,
              'C3prime': C3prime, 's': s, 'y':y}

        end = time.time()
        rt = end - start
        return ct, m, ytest, rt

    def dec_aibbme(self, pp, idstar, ct, dk):
        start = time.time()

        Vidori = Vidj_func(dk, ct, idstar)
        Vid = group.init(ZR, int(str(Vidori)[:2]))
        d2dec = ddec_func(ct, Vid)

        rtagdec, y = rtagdec_func(ct, dk)

        A = A_func(pp, ct, dk, y, d2dec, idstar)
        B = B_func(ct, dk)

        dec_msg = ((A ** (1 / (rtagdec - ct['ctag']))) * ct['C0']) / B
        print('dec message:', dec_msg)
        
        exit(0)
        end = time.time()
        rt = end - start
        return dec_msg, rt

    def rkgen_aibbme(self, n, pp, msk):
        start = time.time()
        sigma = group.random(ZR)
        sprime = group.random(ZR)
        d0 = (pp['g'] ** sigma) + FHash_function(pp, pp['e_gu'] ** sprime)

        tempn = n
        IDnprime = []
        for i in range(tempn):
            IDnprime.append('hello' + str(i + tempn))

        temp = 1
        temp2 = 1
        for i in range(n):
            temp = temp * (alpha + H.hashToZr(IDnprime[i]))
            temp2 = temp2 * ((alpha + H.hashToZr(IDnprime[i])) / (H.hashToZr(IDnprime[i])))

        d1 = pp['g'] ** (sprime * temp)
        d2 = pp['miu1'] ** (-sprime)
        d2prime = msk['miu'] ** (-sprime)

        DK = {'d0': d0, 'd1': d1, 'd2': d2, 'd2prime': d2prime, 'sigma': sigma, 'IDnprime': IDnprime}

        end = time.time()
        ct = end - start
        return DK, ct


def random_array(n):
    array = []
    for i in range(n):
        array.append(group.random(ZR))
    return array


def r_func(n, b, r1n, r2n):
    r = []
    for i in range(n):
        r.append(r1n[i] + b * r2n[i])
    return r


def R_func(n, g, rn):
    Rn = []
    for i in range(n):
        Rn.append(g ** rn[i])
    return Rn


def H0(X):
    X = bytes([a ^ b for (a, b) in zip(X.encode(), bytes.fromhex(mask))])
    return group.hash(X, G2)


def H1(X):
    X = bytes([a ^ b for (a, b) in zip(X.encode(), bytes.fromhex(mask))])
    return group.hash(X, G1)


def H2(X):
    H = Hash(group)
    return H.hashToZr(X)


def H2toint(X):
    H = Hash(group)
    return int(H.hashToZr(X))


def H3(X):
    H = Hash(group)
    return H.hashToZn(X)


def hr1_func(n, h, r1n):
    hr1n = []
    for i in range(n):
        hr1n.append(h ** r1n[i])
    return hr1n


def hr2_func(n, h, r2n):
    hr2n = []
    for i in range(n):
        hr2n.append(h ** r2n[i])
    return hr2n


def rtagsn_func(n):
    rtagsn = []
    for i in range(n):
        rtagsn.append(group.random(ZR))
    return rtagsn


def dk7_func( pp, rtagsn, id, z, ct):
    dk7n = []
    H2idcut = group.init(ZR, int(str(H2(id))[:2]))

    for i in range(1, len(ct['y'])):
        dk71 = (pp['ht1'] ** (rtagsn[i-1])) * pp['hr1n'][i]
        dk72 = pp['hr1n'][0] ** ((H2idcut) ** i)
        dk7n.append((dk71 / dk72) ** z)

    return dk7n


def dk7test_func(pp, id, ct, rtagdec, dk):
    a1 = 0

    for i in range(1, len(ct['y'])):
        a1 = a1 + (dk['rtagsn'][i - 1]) * ct['y'][i]

    a2 = pp['ht1'] ** a1
    a3 = pp['ht1'] ** rtagdec

    return a2

def dk8_func( pp, rtagsn, id, z, ct):
    dk8n = []
    H2idcut = group.init(ZR, int(str(H2(id))[:2]))

    for i in range(1, len(ct['y'])):
        dk81 = (pp['ht2'] ** (rtagsn[i-1])) * pp['hr2n'][i]
        dk82 = pp['hr2n'][0] ** ((H2idcut) ** i)
        dk8n.append((dk81 / dk82) ** z)

    return dk8n

# def dk8_func(pp, rtagsn, id, z, ct):
#     dk8n = []
#     H2idcut = group.init(ZR, int(str(H2(id))[:2]))
#
#     for i in range(1, len(ct['y'])):
#         dk81 = (pp['ht2'] ** (rtagsn[i-1])) * pp['hr2n'][i]
#         dk82 = pp['hr2n'][0] ** ((H2idcut) ** i)
#         dk8n.append((dk81 / dk82) ** z)
#     return dk8n


def id_generate_func(len):
    res = ''.join(random.choices(string.ascii_uppercase +
                                 string.digits, k=len))
    return res


def idj_func(id_num):
    res = []
    for i in range(id_num):
        res.append(id_generate_func(id_len))
    return res


def H2idj_func(idj):
    res = []
    for i in range(len(idj)):
        res.append(H2toint(idj[i]))
    return res


def cut_func(H2idjtemp, k):
    res = []
    for i in range(len(H2idjtemp)):
        res.append(int(str(H2idjtemp[i])[:k]))
    return res


def coeff_func(H2idj):
    res = np.array([1, -H2idj[0]])
    for i in range(len(H2idj) - 1):
        a2 = np.array([1, -H2idj[i + 1]])
        res = np.convolve(res, a2)
    rev = []
    for i in range(len(res)):
        rev.append(res[len(res) - 1 - i ])
    return rev


def bk_func(vidtemp1, d2):
    vidtemp1[0] = vidtemp1[0] + d2
    return vidtemp1


def Vidj_func(dk, ct, idstar):
    a1 = pair(dk['dk3'], ct['C2'])
    a2 = pair(dk['dk2'], H1(idstar))
    a3 = pair(dk['dk1'], ct['C4'])
    res = H3(a1 * a2 * a3)
    return res


# def d2dec_func(ct, Vid):
#     res = 1
#     for i in range(len(ct['idj'])-1):
#         res = res * (ct['bk'][i] * (int(Vid) ** i))
#
#     res = res + int(Vid) ** (len(ct['idj']))
#     return res

def ddec_func(ct, Vid):
    dd = 0
    for i in range(len(ct['bk'])):
        dd = dd + ct['bk'][i] * (Vid ** i)
    res = int2elment(dd)
    return res


# def ddec_func(an, Uid):
#     dd = 0
#     for i in range(len(an)):
#         dd = dd + an[i] * (Uid ** (len(an) - i - 1))
#     res = int2elment(dd)
#     return res

def int2elment(val):
    return group.init(ZR, int(val))


def rtagdec_func(ct, dk):
    H2idjtemp = H2idj_func(ct['idj'])
    H2idj = cut_func(H2idjtemp, 2)
    # print("dec H2idj: ", H2idj)

    ytemp = coeff_func(H2idj)
    y = yelement_func(ytemp)
    print("dec y: ", y)

    res = 0
    for i in range(1, len(y)):
        res = res + (dk['rtagsn'][i-1] * y[i])
    return res, y




def A_func(pp, ct, dk, y, d2dec, idstar):
    A12 = 1
    for i in range(len(y)-1):
        A12 = A12 * (dk['dk7n'][i] ** y[i+1])
    # A12test = A12test_func(pp, ct, dk)
    dk7ny = dk7ny_func(pp, y, dk)
    # test111 =  test111_func(pp, ct, dk, y, d2dec, idstar)

    a1 = pair(ct['C1'], A12)
    A22 = 1
    for i in range(len(y)-1):
        A22 = A22 * (dk['dk8n'][i] ** y[i+1])
    # A22test = A22test_func(pp, ct, dk)

    a2 = pair(ct['C2'], A22)
    a3 = pair(ct['C3'] ** (1 / d2dec), dk['dk6'])
    return a1 * a2 / a3


def dk7ny_func(pp, y, dk):
    a1 = 0
    for i in range(1, len(y)):
        a1 = a1 + (dk['rtagsn'][i-1] * y[i])
    a2 = pp['ht1'] ** a1
    a3 = 0
    for i in range(len(y)):
        a3 = a3 + pp['r1ntest'][i] ** y[i]
    a4 = pp['h'] ** a3
    a5 = (a2 * a4) ** dk["z"]
    return a5

def test111_func(pp, ct, dk, y, d2dec, idstar):
    a1 = 0
    for i in range(1, len(y)):
        a1  =  a1 + pp['r1ntest'][i] * y[i]
    a2 = pp["h"] ** a1
    a3 = a2 /(pp['h'] ** (- pp['r1ntest'][0] * y[0]))

    a4 = 0
    for i in range(len(y)):
        a4 = a4 + pp['r1ntest'][i] * y[i]
    a5 = pp["h"] ** a4

    return a3



def A12test_func(pp, ct, dk):
    a1 = 0
    for i in range(1, len(ct['y'])):
        a1 = a1 + (dk['rtagsn'][i-1] * ct['y'][i])
    a2 = pp['ht1'] ** a1
    a3 = 0
    for i in range(len(ct['y'])):
        a3 = a3 + pp['r1ntest'][i] ** ct['y'][i]
    a4 = pp['h'] ** a3
    a5 = (a2 * a4) ** dk['z']
    return a5

def A22test_func(pp, ct, dk):
    a1 = 0
    for i in range(1, len(ct['y'])):
        a1 = a1 + (dk['rtagsn'][i - 1] * ct['y'][i])
    a2 = pp['ht2'] ** a1
    a3 = 0
    for i in range(len(ct['y'])):
        a3 = a3 + pp['r2ntest'][i] ** ct['y'][i]
    a4 = pp['h'] ** a3
    a5 = (a2 * a4) ** dk['z']
    return a5



def B_func(ct, dk):
    res = pair(ct['C1'], dk['dk4']) * pair(ct['C2'], dk['dk5'])
    return res


def Atest_func(pp, idstar, ct, dk):
    a1 = pair(pp['T'], pp['h'])
    a2 = 0
    for i in range(1, len(ct['y'])):
        a2 = a2 + (dk['rtagsn'][i] * ct['y'][i])
    ktag = a2
    a3 = (ktag - ct['ctag']) * ct['s'] * dk['z']
    print("ktag - ct['ctag']: ", ktag - ct['ctag'])
    return a1 ** a3





def Btest_func(pp, idstar, ct, dk):
    a1 = pp['X']  ** ct["s"]
    a2 = pair(pp['T'], pp['h']) ** (dk['z'] * ct['s'])
    return a1 * a2


def ypositive_func(ytemp):
    res = []
    for i in range(len(ytemp)):
        if ytemp[i] >= 0:
            res.append(group.init(ZR, int(ytemp[i])))
        else:
            res.append(group.init(ZR, int(-ytemp[i])))
    return res

def yelement_func(ytemp):
    res = []
    for i in range(len(ytemp)):
        res.append(group.init(ZR, int(ytemp[i])))
    return res



def C3_func(pp, d2, s, ctag, y):
    res = 1
    for i in range(len(y)):
        res = res * (pp['R'][i] ** y[i])
    res = ((pp["T"] ** ctag) * res) ** (d2 * s)
    return res


def C3prime_func(pp, s, ctag, y):
    res = 1
    for i in range(len(y)):
        res = res * (pp['R'][i] ** y[i])

    res = ((pp["T"] ** ctag) * res) ** s
    return res


def Vid_func(pp, ek, idj, s):
    res = []
    for i in range(len(idj)):
        a1 = H0(idj[i])
        a2 = ek['ekid'] * (pp["gb"] ** s) * (pp['vrho'] ** s)
        res.append(int(H3(pair(a1, a2))))
    return res


def gN_function(n, g, alpha):
    res = []
    for i in range(n):
        res.append(g ** (alpha ** i))
    return res


def generate_random_str(length):
    random_str = ''
    base_str = 'helloworlddfafj23i4jri3jirj23idaf2485644f5551jeri23jeri23ji23'
    for i in range(length):
        random_str += base_str[random.randint(0, length - 1)]
    return random_str


def FHash_function(pp, inputGT):
    equ1 = H.hashToZn(inputGT)
    equ2 = H.hashToZr(equ1)
    res = pp['g'] ** equ2
    return res


def temp_function(Snew, alpha):
    temp1 = 1
    temp2 = 1
    for i in range(len(Snew)):
        temp1 = temp1 * (alpha + H.hashToZr(Snew[i]))
        temp2 = temp2 * (H.hashToZr(Snew[i]))

    return temp1, temp2


def main():
    groupObj = PairingGroup('SS512')
    # n_array = np.arange(5, 30, 5)
    n_array = [10]
    output_txt = './30_aibbme.txt'
    ahnipe = MJ18(groupObj)

    with open(output_txt, 'w+', encoding='utf-8') as f:
        f.write(
            "Seq SetupAveTime       KeygenAveTime      EncAveTime         Dec1AVeTime        RekeygenAveTime    ReencAveTime       Dec2AveTime   " + '\n')

        for i in range(len(n_array)):
            seq = 5
            sttot, kgtot, enctot, dec1tot, rktot, retot, dec2tot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            for j in range(seq):
                n = n_array[i]
                idstar = 'hello123'

                pp, msk, setuptime = ahnipe.setup_aibbme()
                ek, ekgentime = ahnipe.ekgen_aibbme(msk, idstar)
                ct, m, ytest, enctime = ahnipe.enc_aibbme(pp, ek)
                dk, dkgentime = ahnipe.dkgen_aibbme(pp, msk, ct['idj'][0], ct)
                rec_msg, dectime = ahnipe.dec_aibbme(pp, idstar, ct, dk)

                print('\nn, seq:   ', n, j)
                print("m:        ", m)
                print("rec_msg1: ", rec_msg)
                # print("rec_msg2: ", rec_msg2)

                # m_inputkey = group.serialize(m).decode("utf-8")
                # m_outputkey = group.serialize(rec_msg1).decode("utf-8")
                # encrypt(m_inputkey)
                # decrypt(m_outputkey)
                # image.encrypt(m_inputkey)
                # image.decrypt(m_outputkey)

                # sttot, kgtot, enctot, dec1tot, rktot, retot, dec2tot = sttot + setuptime, kgtot + keygen1time + keygen2time, enctot + enctime, dec1tot + dec1time, rktot + rkgentime, retot + reenctime, dec2tot + dec2time

            out0 = str(n).zfill(2)
            out1 = str(format(sttot / float(seq), '.16f'))
            out2 = str(format(kgtot / float(seq), '.16f'))
            out3 = str(format(enctot / float(seq), '.16f'))
            out4 = str(format(dec1tot / float(seq), '.16f'))
            out5 = str(format(rktot / float(seq), '.16f'))
            out6 = str(format(retot / float(seq), '.16f'))
            out7 = str(format(dec2tot / float(seq), '.16f'))
            f.write(out0 + '  ' + out1 + ' ' + out2 + ' ' + out3 + ' ' + out4 + ' ' + out5 + ' ' + out6 + ' ' + out7)
            f.write('\n')


if __name__ == "__main__":
    main()

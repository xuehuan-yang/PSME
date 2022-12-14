"""
# H0: {0, 1}* -> G1   H_prime()
# H1: {0, 1}* -> G0   H_prime()

# from charm.toolbox.hash_module import Hash
# H = Hash(group)
# H.hashToZr()
# H.hashtoZn()
# H2: {0 ,1}* -> Zp  H.hashToZr()
# H3: GT -> Zp       H.hashtoZn()

"""
from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.ABEncMultiAuth import ABEncMultiAuth
import time
import numpy as np
from charm.toolbox.hash_module import Hash
import sys

sys.path.append('../')
from common.image import *
from common.msp import *
from common.id import *
from common.parse import *


class MJ18(ABEncMultiAuth):
    def __init__(self, groupObj, verbose=False):
        ABEncMultiAuth.__init__(self)
        global group, ahnipe, util, H, mask, l
        group = groupObj
        util = MSP(group, verbose=False)
        H = Hash(group)
        mask = 'ed27dbfb02752e0e16bc4502d6c732bc5f1cc92ba19b2d93a4e95c597ca42753e93550b52f82b6c13fb8cc0c2fc64487'
        l = 100

    def setup_aibbme(self):
        start = time.time()

        g, v = group.random(G1), group.random(G1)
        h = group.random(G2)
        r1n, r2n = random_array(l + 1), random_array(l + 1)
        t1, t2, beta1, beta2, alpha, rho = group.random(ZR, 6)
        b, tau = group.random(ZR), group.random(ZR)

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

    def enc_aibbme(self, pp, ek, n):
        start = time.time()

        idj = idj_func(n)
        H2idj = H2idj_func(idj)

        y_ = coeff_func(H2idj)
        y = yelement_func(y_)

        m = group.random(GT)
        s, d2, ctag = group.random(ZR, 3)
        print("enc d2:   ", d2)

        C0 = m * pp['X'] ** s
        C1 = pp["g"] ** s
        C2 = pp["gb"] ** s
        C3 = C3_func(pp, d2, s, ctag, y)
        C4 = pp['v'] ** s
        Vid = Vid_func(pp, ek, idj, s)

        bk = bk_func(coeff_func(Vid), d2)

        ct = {'C0': C0, 'C1': C1, 'C2': C2, 'C3': C3, 'C4': C4, 'ctag': ctag, 'bk': bk, 'idj': idj, 'd2': d2,
              's': s, 'y': y}

        end = time.time()
        rt = end - start
        return ct, m, rt

    def dec_aibbme(self, idstar, ct, dk):
        start = time.time()

        Vid = Vidj_func(dk, ct, idstar)
        d2dec = ddec_func(ct, Vid)
        print("dec d2:   ", d2dec)
        rtagdec, y = rtagdec_func(ct, dk)
        A = A_func(ct, dk, y, d2dec)
        B = B_func(ct, dk)
        dec_msg = ((A ** (1 / (rtagdec - ct['ctag']))) * ct['C0']) / B

        end = time.time()
        rt = end - start
        return dec_msg, rt


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
    return int2elment(int(H.hashToZr(X)))


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


def dk7_func(pp, rtagsn, id, z, ct):
    dk7n = []
    H2idcut = int2elment(int(H2(id)))

    for i in range(1, len(ct['y'])):
        dk71 = (pp['ht1'] ** (rtagsn[i - 1])) * pp['hr1n'][i]
        dk72 = pp['hr1n'][0] ** (H2idcut ** group.init(ZR, i))
        dk7n.append((dk71 / dk72) ** z)
    return dk7n


def dk8_func(pp, rtagsn, id, z, ct):
    dk8n = []
    H2idcut = int2elment(int(H2(id)))

    for i in range(1, len(ct['y'])):
        dk81 = (pp['ht2'] ** (rtagsn[i - 1])) * pp['hr2n'][i]
        dk82 = pp['hr2n'][0] ** (H2idcut ** group.init(ZR, i))
        dk8n.append((dk81 / dk82) ** z)
    return dk8n


def idj_func(n):
    res = []
    for i in range(n):
        res.append(id_generator(n))
    return res


def H2idj_func(idj):
    res = []
    for i in range(len(idj)):
        res.append(H2toint(idj[i]))
    return res


def coeff_func(H2idj):
    res = np.array([1, -H2idj[0]])
    for i in range(len(H2idj) - 1):
        a2 = np.array([1, -H2idj[i + 1]])
        res = np.convolve(res, a2)
    rev = []
    for i in range(len(res)):
        rev.append(res[len(res) - 1 - i])
    return rev


def bk_func(vid, d2):
    vid[0] = vid[0] + int(d2)
    return vid


def Vidj_func(dk, ct, idstar):
    a1 = pair(dk['dk3'], ct['C2'])
    a2 = pair(dk['dk2'], H1(idstar))
    a3 = pair(dk['dk1'], ct['C4'])
    res = H3(a1 * a2 * a3)
    return int2elment(int(res))


def ddec_func(ct, Vid):
    dd = group.init(ZR, 0)
    for i in range(len(ct['bk'])):
        dd = dd + ct['bk'][i] * (Vid ** group.init(ZR, i))
    res = int2elment(dd)
    return res


def int2elment(val):
    return group.init(ZR, int(val))


def rtagdec_func(ct, dk):
    H2idj = H2idj_func(ct['idj'])

    y_ = coeff_func(H2idj)
    y = yelement_func(y_)
    res = group.init(ZR, 0)
    for i in range(1, len(y)):
        res = res + (dk['rtagsn'][i - 1] * y[i])
    return res, y


def A_func(ct, dk, y, d2dec):
    A12 = group.init(ZR, 1)
    for i in range(len(y) - 1):
        A12 = A12 * (dk['dk7n'][i] ** y[i + 1])

    a1 = pair(ct['C1'], A12)
    A22 = group.init(ZR, 1)
    for i in range(len(y) - 1):
        A22 = A22 * (dk['dk8n'][i] ** y[i + 1])

    a2 = pair(ct['C2'], A22)
    a3 = pair(ct['C3'] ** (1 / d2dec), dk['dk6'])
    return a1 * a2 / a3


def B_func(ct, dk):
    res = pair(ct['C1'], dk['dk4']) * pair(ct['C2'], dk['dk5'])
    return res


def yelement_func(dir):
    res = []
    for i in range(len(dir)):
        res.append(group.init(ZR, int(dir[i])))
    return res


def C3_func(pp, d2, s, ctag, y):
    res = group.init(ZR, 1)
    for i in range(len(y)):
        res = res * (pp['R'][i] ** y[i])
    res = ((pp["T"] ** ctag) * res) ** (d2 * s)
    return res


def Vid_func(pp, ek, idj, s):
    res = []
    for i in range(len(idj)):
        a1 = H0(idj[i])
        a2 = ek['ekid'] * (pp["gb"] ** s) * (pp['vrho'] ** s)
        res.append(int2elment(int(H3(pair(a1, a2)))))
    return res




def main():
    groupObj = PairingGroup('SS512')
    n_array = np.arange(5, 55, 5)
    n_array = np.insert(n_array, 0, 1)
    output_txt= output_func('30_aibbme')
    ahnipe = MJ18(groupObj)

    with open(output_txt, 'w+', encoding='utf-8') as f:
        f.write(
            "Seq SetupAveTime       ekgenAveTime       EncAveTime         dkgenAveTime       decAveTime        " + '\n')

        for i in range(len(n_array)):
            seq = seq_func()
            sttot, ekgentot, enctot, dkgentot, dectot = 0.0, 0.0, 0.0, 0.0, 0.0
            for j in range(seq):
                n = n_array[i]
                idstar = id_generator(n)
                print('\n')
                print('n, seq:   ', n, j)

                pp, msk, setuptime = ahnipe.setup_aibbme()
                ek, ekgentime = ahnipe.ekgen_aibbme(msk, idstar)
                ct, m, enctime = ahnipe.enc_aibbme(pp, ek, n)
                dk, dkgentime = ahnipe.dkgen_aibbme(pp, msk, ct['idj'][0], ct)
                rec_msg, dectime = ahnipe.dec_aibbme(idstar, ct, dk)

                print("m:        ", m)
                print("rec_msg1: ", rec_msg)

                m_inputkey = group.serialize(m).decode("utf-8")
                m_outputkey = group.serialize(rec_msg).decode("utf-8")
                encrypt(m_inputkey)
                decrypt(m_outputkey)

                sttot, ekgentot, enctot, dkgentot, dectot = sttot + setuptime, ekgentot + ekgentime, enctot + enctime, dkgentot + dkgentime, dectot + dectime

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

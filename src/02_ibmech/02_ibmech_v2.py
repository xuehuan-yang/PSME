from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, G2, GT, pair
from charm.toolbox.matrixops import *
from charm.core.crypto.cryptobase import *
from charm.toolbox.IBEnc import IBEnc

debug = True


class IBE_Chen12_z(IBEnc):

    def __init__(self, groupObj):
        IBEnc.__init__(self)
        global group
        group = groupObj

    def setup(self):
        g1 = group.random(G1)
        g2 = group.random(G2)
        alpha = group.random(ZR)
        eta = group.random(ZR)

        d11, d12, d13, d14, d15, d16, d17, d18 = group.random(ZR, 8)
        d21, d22, d23, d24, d25, d26, d27, d28 = group.random(ZR, 8)
        d31, d32, d33, d34, d35, d36, d37, d38 = group.random(ZR, 8)
        d41, d42, d43, d44, d45, d46, d47, d48 = group.random(ZR, 8)
        d51, d52, d53, d54, d55, d56, d57, d58 = group.random(ZR, 8)
        d61, d62, d63, d64, d65, d66, d67, d68 = group.random(ZR, 8)
        d71, d72, d73, d74, d75, d76, d77, d78 = group.random(ZR, 8)
        d81, d82, d83, d84, d85, d86, d87, d88 = group.random(ZR, 8)

        one = group.random(ZR)

        [D11, D12, D13, D14, D15, D16, D17, D18] = GaussEliminationinGroups(
            [[d11, d12, d13, d14, d15, d16, d17, d18, one],
             [d21, d22, d23, d24, d25, d26, d27, d28, group.init(ZR, 0)],
             [d31, d32, d33, d34, d35, d36, d37, d38, group.init(ZR, 0)],
             [d41, d42, d43, d44, d45, d46, d47, d48, group.init(ZR, 0)],
             [d51, d52, d53, d54, d55, d56, d57, d58, group.init(ZR, 0)],
             [d61, d62, d63, d64, d65, d66, d67, d68, group.init(ZR, 0)],
             [d71, d72, d73, d74, d75, d76, d77, d78, group.init(ZR, 0)],
             [d81, d82, d83, d84, d85, d86, d87, d88, group.init(ZR, 0)]])

        [D21, D22, D23, D24, D25, D26, D27, D28] = GaussEliminationinGroups(
            [[d11, d12, d13, d14, d15, d16, d17, d18, group.init(ZR, 0)],
             [d21, d22, d23, d24, d25, d26, d27, d28, one],
             [d31, d32, d33, d34, d35, d36, d37, d38, group.init(ZR, 0)],
             [d41, d42, d43, d44, d45, d46, d47, d48, group.init(ZR, 0)],
             [d51, d52, d53, d54, d55, d56, d57, d58, group.init(ZR, 0)],
             [d61, d62, d63, d64, d65, d66, d67, d68, group.init(ZR, 0)],
             [d71, d72, d73, d74, d75, d76, d77, d78, group.init(ZR, 0)],
             [d81, d82, d83, d84, d85, d86, d87, d88, group.init(ZR, 0)]])

        [D31, D32, D33, D34, D35, D36, D37, D38] = GaussEliminationinGroups(
            [[d11, d12, d13, d14, d15, d16, d17, d18, group.init(ZR, 0)],
             [d21, d22, d23, d24, d25, d26, d27, d28, group.init(ZR, 0)],
             [d31, d32, d33, d34, d35, d36, d37, d38, one],
             [d41, d42, d43, d44, d45, d46, d47, d48, group.init(ZR, 0)],
             [d51, d52, d53, d54, d55, d56, d57, d58, group.init(ZR, 0)],
             [d61, d62, d63, d64, d65, d66, d67, d68, group.init(ZR, 0)],
             [d71, d72, d73, d74, d75, d76, d77, d78, group.init(ZR, 0)],
             [d81, d82, d83, d84, d85, d86, d87, d88, group.init(ZR, 0)]])

        [D41, D42, D43, D44, D45, D46, D47, D48] = GaussEliminationinGroups(
            [[d11, d12, d13, d14, d15, d16, d17, d18, group.init(ZR, 0)],
             [d21, d22, d23, d24, d25, d26, d27, d28, group.init(ZR, 0)],
             [d31, d32, d33, d34, d35, d36, d37, d38, group.init(ZR, 0)],
             [d41, d42, d43, d44, d45, d46, d47, d48, one],
             [d51, d52, d53, d54, d55, d56, d57, d58, group.init(ZR, 0)],
             [d61, d62, d63, d64, d65, d66, d67, d68, group.init(ZR, 0)],
             [d71, d72, d73, d74, d75, d76, d77, d78, group.init(ZR, 0)],
             [d81, d82, d83, d84, d85, d86, d87, d88, group.init(ZR, 0)]])

        # generate public parameters.
        PP2 = (pair(g1, g2)) ** (alpha * one)
        PP3 = (pair(g1, g2)) ** (eta * one)
        gd11, gd12, gd13, gd14, gd15, gd16, gd17, gd18 = g1 ** d11, g1 ** d12, g1 ** d13, g1 ** d14, g1 ** d15, g1 ** d16, g1 ** d17, g1 ** d18
        gd21, gd22, gd23, gd24, gd25, gd26, gd27, gd28 = g1 ** d21, g1 ** d22, g1 ** d23, g1 ** d24, g1 ** d25, g1 ** d26, g1 ** d27, g1 ** d28
        gd31, gd32, gd33, gd34, gd35, gd36, gd37, gd38 = g1 ** d31, g1 ** d32, g1 ** d33, g1 ** d34, g1 ** d35, g1 ** d36, g1 ** d37, g1 ** d38
        gd41, gd42, gd43, gd44, gd45, gd46, gd47, gd48 = g1 ** d41, g1 ** d42, g1 ** d43, g1 ** d44, g1 ** d45, g1 ** d46, g1 ** d47, g1 ** d48
        # gd51, gd52, gd53, gd54, gd55, gd56, gd57, gd58 = g1 ** d51, g1 ** d52, g1 ** d53, g1 ** d54, g1 ** d55, g1 ** d56, g1 ** d57, g1 ** d58
        # gd61, gd62, gd63, gd64, gd65, gd66, gd67, gd68 = g1 ** d61, g1 ** d62, g1 ** d63, g1 ** d64, g1 ** d65, g1 ** d66, g1 ** d67, g1 ** d68
        # gd71, gd72, gd73, gd74, gd75, gd76, gd77, gd78 = g1 ** d71, g1 ** d72, g1 ** d73, g1 ** d74, g1 ** d75, g1 ** d76, g1 ** d77, g1 ** d78
        # gd81, gd82, gd83, gd84, gd85, gd86, gd87, gd88 = g1 ** d81, g1 ** d82, g1 ** d83, g1 ** d84, g1 ** d85, g1 ** d86, g1 ** d87, g1 ** d88

        pk = {'PP2': PP2, 'PP3': PP3,
              'gd11': gd11, 'gd12': gd12, 'gd13': gd13, 'gd14': gd14, 'gd15': gd15, 'gd16': gd16, 'gd17': gd17,'gd18': gd18,
              'gd21': gd21, 'gd22': gd22, 'gd23': gd23, 'gd24': gd24, 'gd25': gd25, 'gd26': gd26, 'gd27': gd27,'gd28': gd28, }

        msk = {'alpha': alpha, 'eta': eta, 'g1': g1, 'g2': g2,
               'gd31': gd31, 'gd32': gd32, 'gd33': gd33, 'gd34': gd34, 'gd35': gd35, 'gd36': gd36, 'gd37': gd37, 'gd38': gd38,
               'gd41': gd41, 'gd42': gd42, 'gd43': gd43, 'gd44': gd44, 'gd45': gd45, 'gd46': gd46, 'gd47': gd47, 'gd48': gd48,
               'D11': D11, 'D12': D12, 'D13': D13, 'D14': D14, 'D15': D15, 'D16': D16, 'D17': D17, 'D18': D18,
               'D21': D21, 'D22': D22, 'D23': D23, 'D24': D24, 'D25': D25, 'D26': D26, 'D27': D27, 'D28': D28,
               'D31': D31, 'D32': D32, 'D33': D33, 'D34': D34, 'D35': D35, 'D36': D36, 'D37': D37, 'D38': D38,
               'D41': D41, 'D42': D42, 'D43': D43, 'D44': D44, 'D45': D45, 'D46': D46, 'D47': D47, 'D48': D48,
               }
        if (debug):
            # print("Public parameters...")
            group.debug(pk)
            # print("Secret parameters...")
            group.debug(msk)
        return (pk, msk)

    def skgen(self, msk, ID):
        _ID = group.hash(ID)
        r = group.random(ZR)
        ek_id1 = msk['gd31'] ** (msk['eta'] + r * _ID) / (msk['gd41'] ** r)
        ek_id2 = msk['gd32'] ** (msk['eta'] + r * _ID) / (msk['gd42'] ** r)
        ek_id3 = msk['gd33'] ** (msk['eta'] + r * _ID) / (msk['gd43'] ** r)
        ek_id4 = msk['gd34'] ** (msk['eta'] + r * _ID) / (msk['gd44'] ** r)
        ek_id5 = msk['gd35'] ** (msk['eta'] + r * _ID) / (msk['gd45'] ** r)
        ek_id6 = msk['gd36'] ** (msk['eta'] + r * _ID) / (msk['gd46'] ** r)
        ek_id7 = msk['gd37'] ** (msk['eta'] + r * _ID) / (msk['gd47'] ** r)
        ek_id8 = msk['gd38'] ** (msk['eta'] + r * _ID) / (msk['gd48'] ** r)

        k = {'ek_id1': ek_id1, 'ek_id2': ek_id2, 'ek_id3': ek_id3,
             'ek_id4': ek_id4, 'ek_id5': ek_id5, 'ek_id6': ek_id6, 'ek_id7': ek_id7,
             'ek_id8': ek_id8, }

        if (debug):
            # print("Generate User SK...")
            group.debug(k)
        return k

    def rkgen(self, pk, msk, IDstar):

        _IDstar = group.hash(IDstar)
        s, s1, s2 = group.random(ZR, 3)

        sk_id11 = msk['g2'] ** (msk['D11'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D21'] + s * msk['D31'])
        sk_id12 = msk['g2'] ** (msk['D12'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D22'] + s * msk['D32'])
        sk_id13 = msk['g2'] ** (msk['D13'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D23'] + s * msk['D33'])
        sk_id14 = msk['g2'] ** (msk['D14'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D24'] + s * msk['D34'])
        sk_id15 = msk['g2'] ** (msk['D15'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D25'] + s * msk['D35'])
        sk_id16 = msk['g2'] ** (msk['D16'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D26'] + s * msk['D36'])
        sk_id17 = msk['g2'] ** (msk['D17'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D27'] + s * msk['D37'])
        sk_id18 = msk['g2'] ** (msk['D18'] * (msk['alpha'] + s1 * _IDstar) - s1 * msk['D28'] + s * msk['D38'])

        sk_id21 = msk['g2'] ** (s2 * (_IDstar * msk['D11'] - msk['D21']) + s * msk['D41'])
        sk_id22 = msk['g2'] ** (s2 * (_IDstar * msk['D12'] - msk['D22']) + s * msk['D42'])
        sk_id23 = msk['g2'] ** (s2 * (_IDstar * msk['D13'] - msk['D23']) + s * msk['D43'])
        sk_id24 = msk['g2'] ** (s2 * (_IDstar * msk['D14'] - msk['D24']) + s * msk['D44'])
        sk_id25 = msk['g2'] ** (s2 * (_IDstar * msk['D15'] - msk['D25']) + s * msk['D45'])
        sk_id26 = msk['g2'] ** (s2 * (_IDstar * msk['D16'] - msk['D26']) + s * msk['D46'])
        sk_id27 = msk['g2'] ** (s2 * (_IDstar * msk['D17'] - msk['D27']) + s * msk['D47'])
        sk_id28 = msk['g2'] ** (s2 * (_IDstar * msk['D18'] - msk['D28']) + s * msk['D48'])

        sk_id3 =  pk['PP3'] ** s

        sk = {'sk_id11': sk_id11, 'sk_id12': sk_id12, 'sk_id13': sk_id13,'sk_id14': sk_id14,
             'sk_id15': sk_id15, 'sk_id16': sk_id16, 'sk_id17': sk_id17, 'sk_id18': sk_id18,
             'sk_id21': sk_id21, 'sk_id22': sk_id22, 'sk_id23': sk_id23, 'sk_id24': sk_id24,
             'sk_id25': sk_id25, 'sk_id26': sk_id26, 'sk_id27': sk_id27, 'sk_id28': sk_id28,
             "sk_id3":sk_id3, 'stest': s }

        if (debug):
            # print("Generate User SK...")
            group.debug(sk)
        return sk

    def encrypt(self, pk, ek, IDstar, M):
        z = group.random(ZR)
        _IDstar = group.hash(IDstar)
        C0 = (pk['PP2'] ** z) * M

        C11 = (pk['gd11'] ** z) * (pk['gd21'] ** (z * _IDstar)) * (ek['ek_id1'])
        C12 = (pk['gd12'] ** z) * (pk['gd22'] ** (z * _IDstar)) * (ek['ek_id2'])
        C13 = (pk['gd13'] ** z) * (pk['gd23'] ** (z * _IDstar)) * (ek['ek_id3'])
        C14 = (pk['gd14'] ** z) * (pk['gd24'] ** (z * _IDstar)) * (ek['ek_id4'])
        C15 = (pk['gd15'] ** z) * (pk['gd25'] ** (z * _IDstar)) * (ek['ek_id5'])
        C16 = (pk['gd16'] ** z) * (pk['gd26'] ** (z * _IDstar)) * (ek['ek_id6'])
        C17 = (pk['gd17'] ** z) * (pk['gd27'] ** (z * _IDstar)) * (ek['ek_id7'])
        C18 = (pk['gd18'] ** z) * (pk['gd28'] ** (z * _IDstar)) * (ek['ek_id8'])

        ct = {'C0': C0, 'C11': C11, 'C12': C12, 'C13': C13, 'C14': C14,
              'C15': C15, 'C16': C16, 'C17': C17, 'C18': C18, 'ztest': z}

        if (debug):
            # print('\nEncrypt...')
            group.debug(ct)
        return ct

    def decrypt(self, pk, sk,  ID,  ct):
        _ID = group.hash(ID)

        a1 = (  pair(ct['C11'], sk['sk_id11'] * (sk['sk_id21'] ** _ID))
              * pair(ct['C12'], sk['sk_id12'] * (sk['sk_id22'] ** _ID))
              * pair(ct['C13'], sk['sk_id13'] * (sk['sk_id23'] ** _ID))
              * pair(ct['C14'], sk['sk_id14'] * (sk['sk_id24'] ** _ID))
              * pair(ct['C15'], sk['sk_id15'] * (sk['sk_id25'] ** _ID))
              * pair(ct['C16'], sk['sk_id16'] * (sk['sk_id26'] ** _ID))
              * pair(ct['C17'], sk['sk_id17'] * (sk['sk_id27'] ** _ID))
              * pair(ct['C18'], sk['sk_id18'] * (sk['sk_id28'] ** _ID))
              )

        # a1test = (pk['PP3'] ** (sk['stest'])) * (pk['PP2'] ** ct['ztest'])
        Mprime = ct['C0'] * sk['sk_id3'] / a1

        return Mprime


def main():
    group = PairingGroup('SS512', secparam=1024)
    ibe = IBE_Chen12_z(group)
    (pk, msk) = ibe.setup()
    ID = 'user@email.com'
    IDstar = 'name@gmail.com'

    ek = ibe.skgen(msk, ID)
    sk = ibe.rkgen(pk, msk, IDstar)
    msg = group.random(GT)
    print("ori message: ", msg)
    ct = ibe.encrypt(pk, ek, IDstar, msg)
    decryptedMSG = ibe.decrypt(pk, sk, ID,  ct)
    print("dec message: ", decryptedMSG)
    print(decryptedMSG == msg)


if __name__ == '__main__':
    debug = True
    main()
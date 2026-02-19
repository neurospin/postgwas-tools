#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import scipy.stats as stats
import argparse

def main():
    parser = argparse.ArgumentParser(
        description=("Extract mostest results.\n"))
    parser.add_argument("--bim", type=str,
                        help="Path to .bim file (reference set of SNPs)")
    parser.add_argument("--fname", type=str,
                        help="Prefix of .mat files output by mostest.m")
    parser.add_argument("--out", default=None,
                        help="Optional suffix for output files (defaults to fname)")

    # Parse CLI arguments
    args = parser.parse_args()

    # Apply logic: use fname as out if not provided
    bim_file = args.bim
    fname = args.fname
    out = args.out if args.out else fname

    # read .bim file (reference set of SNPs)
    print('Load {}...'.format(bim_file))
    bim = pd.read_csv(bim_file, sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
    del bim['GP']

    mat = sio.loadmat(fname + '.mat')

    """# save figures with  QQ plots MOST and minP
    print('Generate {}.plot.png...'.format(out))
    with np.errstate(divide='ignore'):
        df = pd.DataFrame({
        'minp_x':np.transpose(mat['hv_maxlogpvecs'].flatten()),
        'minp_y1':np.transpose(-np.log10(1-mat['chc_maxlogpvecs'])).flatten(),
        'minp_y2':np.transpose(-np.log10(1-mat['cdf_minpvecs'])).flatten(),
        'most_x':mat['hv_mostvecs'].flatten(),
        'most_y1':np.transpose(-np.log10(1-mat['chc_mostvecs'])).flatten(),
        'most_y2':np.transpose(-np.log10(1-mat['cdf_mostvecs'])).flatten() })
    df.to_csv(out + '.plot.csv',index=False, sep='\t')
    plt.figure(figsize=(20, 10), dpi=100)
    plt.subplot(2,4,1)
    plt.plot(df['minp_x'], df['minp_y1'])
    plt.plot(df['minp_x'], df['minp_y2'])
    plt.legend(['data (null)', 'minP (model)'])
    plt.title('minP')
    plt.subplot(2,4,2)
    plt.plot(df['most_x'], df['most_y1'])
    plt.plot(df['most_x'], df['most_y2'])
    plt.title('MOSTest')
    plt.legend(['data (null)', 'MOSTest (model)'])
    plt.savefig(out + '.plot.png', bbox_inches='tight')"""

    # generate .sumstats files, compatible with FUMA
    print('Generate {}.***.sumstats files...'.format(out))
    for test in ['most_log10pval_orig', 'minp_log10pval_orig']:# 'minp_log10pval_orig','minp_log10pval_perm', 'most_log10pval_perm'
        bim['PVAL'] = np.power(10, -mat[test].flatten())
        bim['Z_FAKE'] = -stats.norm.ppf(bim['PVAL'].values*0.5) #*effect_sign.astype(np.float64) - effect size not available from MOSTest and minP
        bim['N'] = mat['nvec']
        bim.to_csv('{}.{}.sumstats'.format(out, test.replace('_log10pval', '')), sep='\t', na_rep="NA", index=False)

    print('Done.')

if __name__ == '__main__':
    main()

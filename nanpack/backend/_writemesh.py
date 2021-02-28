
def _save_metrics_2d(x, y, xix, xiy, etax, etay, jj, f_name):
    iM, jM = x.shape
    outf = open(f_name, "w")
    print(f"IMAX= {iM}, JMAX= {jM}", file=outf)
    print(f'{"X":^11} {"Y":^11} {"XiX":^11} {"XiY":^11} {"EtaX":^11}\
 {"EtaY":^11} {"JJ":^11}', file=outf)
    for i in range(0, iM):
        for j in range(0, jM):
            print(f"{x[i][j]:>11.6f} {y[i][j]:>11.6f} {xix[i][j]:>11.6f}\
 {xiy[i][j]:>11.6f} {etax[i][j]:>11.6f} {etay[i][j]:>11.6f}\
 {jj[i][j]:>11.6f}", file=outf)

    outf.close()


def _save_metrics_1d(x, xix, jj, f_name):
    iM, = x.shape
    outf = open(f_name, "w")
    print(f"IMAX= {iM}", file=outf)
    print(f'{"X":^11} {"XiX":^11} {"JJ":^11}', file=outf)
    for i in range(0, iM):
        print(f"{x[i]:>11.6f} {xix[i]:>11.6f} {jj[i]:>11.6f}",
              file=outf)

    outf.close()

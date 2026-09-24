from run_case import *
for D in [200,500,2000]:
    ic, peaks, fits, res, truth = run(D, 0)
    nrl=sum(p.resolution_limited for p in peaks)
    print(f"\n=== D_true={D} nm: {len(peaks)} peaks, {nrl} resolution-limited; ratios:", [round(p.broadening_ratio,2) for p in peaks])
    print(ve.ResultsFormatter.format_comparison(res).split('--- Interp')[0])
    if D==2000:
        print(ve.ResultsFormatter.format_scherrer(res['Scherrer']).split('--- Method')[0][-600:])

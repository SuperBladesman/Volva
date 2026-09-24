import numpy as np, volva_engine as ve, warnings, logging, time, os
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
tt=np.arange(10,80,0.02); y=100+1000*np.exp(-((tt-30)/0.1)**2)
def w(name, text): open(name,'w').write(text); r=ve.XRDLoader.load(name); print(f"{name:28s} -> {r[2]}")
body="\n".join(f"{a:.4f} {b:.1f}" for a,b in zip(tt,y))
w("hdr_plain.xy","2theta intensity\n"+body)
w("hdr_hash.xy","# 2theta intensity\n"+body)
w("csv_hdr.csv","2Theta,Counts\n"+"\n".join(f"{a:.4f},{b:.1f}" for a,b in zip(tt,y)))
w("euro_decimal.txt","\n".join(f"{a:.4f};{b:.1f}".replace('.',',') for a,b in zip(tt,y)))
w("three_col.xy","\n".join(f"{a:.4f} {b:.1f} {np.sqrt(b):.2f}" for a,b in zip(tt,y)))
w("trailing_text.xy",body+"\nEND\n")
w("bruker_raw_like.txt","[Measurement]\nWavelength=1.5406\n[Data]\n"+body)
# real sample file supplied with Volva
t,i,m=ve.XRDLoader.load("sample_xrd_data.txt"); print(m, t.min(), t.max(), np.median(np.diff(t)))
# SNIP speed
big=np.random.poisson(200,size=11500).astype(float)
t0=time.time(); ve.BackgroundProcessor.snip(big,40); print("SNIP 11500 pts x40 it:", round(time.time()-t0,2),"s")
t0=time.time(); ve.BackgroundProcessor.snip(big,150); print("SNIP 11500 pts x150 it:", round(time.time()-t0,2),"s")

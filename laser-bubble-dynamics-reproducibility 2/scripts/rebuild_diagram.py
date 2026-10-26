from pathlib import Path
import html,io
from PIL import ImageFont
import fitz
R=Path(__file__).resolve().parents[1];O=R/'outputs/diagram_rebuilt';O.mkdir(parents=True,exist_ok=True);W,H=2000,1630
from matplotlib import font_manager
freg=font_manager.findfont('DejaVu Sans');fbold=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="146.7mm" viewBox="0 0 {W} {H}">', '<title>Laser-driven bubble model implementations and proposed extensions</title>', '<desc>Three separate implemented source-model variants, with diagnostic-only optical fractions distinguished from feedback and proposed physical extensions shown with dashed outlines.</desc>', '<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L10,5 L0,10 z" fill="#536578"/></marker></defs>',f'<rect width="{W}" height="{H}" fill="white"/>']
def rect(x,y,w,h,fill,stroke,sw=2,dash=None,rx=12):
 if dash:
  parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}"/>')
  for yy in [y,y+h]:
   for xx in range(int(x),int(x+w),17):parts.append(f'<line x1="{xx}" y1="{yy}" x2="{min(xx+10,x+w)}" y2="{yy}" stroke="{stroke}" stroke-width="{sw}"/>')
  for xx in [x,x+w]:
   for yy in range(int(y),int(y+h),17):parts.append(f'<line x1="{xx}" y1="{yy}" x2="{xx}" y2="{min(yy+10,y+h)}" stroke="{stroke}" stroke-width="{sw}"/>')
 else:
  parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def text(x,y,s,size=28,bold=False,fill='#172B3D',anchor='start'):
 parts.append(f'<text x="{x}" y="{y}" font-family="Arial" font-size="{size}" font-weight="'+('bold' if bold else 'normal')+f'" fill="{fill}" text-anchor="{anchor}">{html.escape(s)}</text>')
def line(x1,y1,x2,y2,stroke='#536578',width=3,arrow=False):
 parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2-9 if arrow else x2}" y2="{y2}" stroke="{stroke}" stroke-width="{width}"/>')
 if arrow:parts.append(f'<polygon points="{x2},{y2} {x2-11},{y2-7} {x2-11},{y2+7}" fill="{stroke}"/>')
def card(x,y,title,lines,col,bg):
 rect(x,y,444,222,bg,col)
 font=ImageFont.truetype(fbold,30)
 assert font.getlength(title)<400,(title,font.getlength(title))
 text(x+22,y+47,title,30,True,col)
 for j,s in enumerate(lines):
  assert ImageFont.truetype(freg,28).getlength(s)<402,(s,ImageFont.truetype(freg,28).getlength(s))
  text(x+22,y+96+j*40,s,28)
def lane(y,letter,title,col,bg,content):
 rect(44,y,1912,53,bg,col,1.4,rx=8);text(65,y+37,letter+'   '+title,32,True,col)
 cy=y+73
 for i,(h,ls) in enumerate(content):
  x=44+489*i;card(x,cy,h,ls,col,bg)
  if i<3:line(x+452,cy+111,x+480,cy+111,arrow=True)
# Context and status key
text(44,43,'Experimental motivation',31,True)
text(44,85,'1064 nm irradiation of Fe₃O₄ beads in liquid',29)
text(44,124,'Observed dark-region growth, contraction and episodic optical emission',28)
line(1470,41,1537,41,'#536578',3);text(1555,50,'Implemented proxy',25)
[line(xx,84,min(xx+10,1537),84,'#75818C',3) for xx in range(1470,1537,17)];text(1555,93,'Further physics / tests',25)
lane(162,'A','Single-bubble source models (Figures 2–3)','#156181','#EFF7FA',[
 ('Laser gate and energy',['60 s ON / 30 s OFF','Stored-energy proxy E','Fig. 3 loading affects source']),
 ('Radial calculation',['Spring force and damping','Energy-to-pressure mapping','Radius and velocity clipped']),
 ('Event rule',['R < R_c and v < 0','Armed + 0.5 s cooldown','Rearm when R > 1.2R_c']),
 ('Flash and diagnostics',['Assigned event energy','Exponential flash response','Fig. 3: T_diag and A_diag'])])
rect(44,477,1912,72,'#F8FBFD','#B4CEDB',1.3,rx=7)
text(66,507,'Figure 3 optics: τ = τ₀(1 + 4M_np),  T_diag = exp(−τ),  A_diag = 1 − T_diag',27)
text(66,538,'These absorption-only fractions are diagnostics and do not enter the energy update.',27,True,'#156181')
lane(577,'B','Two-reservoir ensemble with logistic activation (Figure 5)','#8C5A14','#FFFAED',[
 ('Absorbed laser power',['P_abs = L P_in f_abs(R)','Effective local power scale','80 µs gate, 75% ON']),
 ('Two thermal reservoirs',['Separate T_e, T_i balances','Heating, exchange, cooling','Logistic activation x_e(T_e)']),
 ('Pressure and events*',['Ideal-mixture drive pressure','Instantaneous coupling','Downward radius crossing']),
 ('Camera response',['Burst photon-source proxy','20 µs detector response','Saturation in arbitrary units'])])
text(66,918,'T_e is a modeled reservoir temperature. T_i evolves internally. The source does not calculate a spectrum.',27,fill='#755121')
lane(952,'C','Extended ensemble with Saha and loading feedback (Figure 6)','#63529A','#F5F2FB',[
 ('Loading feedback',['Dimensionless M in [0, 10]','Absorption factor 1 + 0.25M','Pressure-driven loading law']),
 ('Thermal and Saha',['Separate T_e, T_i balances','Effective Saha x_ion(n,T_e)','Effective E_ion = 12.6 eV']),
 ('Pressure interactions*',['Instantaneous coupling','t_arr = t_c + r/c_l','c_l = 1500 m/s; decay 3 µs']),
 ('Events and detector',['Radius-crossing trigger','Emission factor 1 + 0.35⟨M⟩','Camera in arbitrary units'])])
text(66,1297,'Figures 3 and 6 include pressure-impulse or threshold damage indicators, without resolved membrane mechanics.',26,fill='#5F527A')
text(44,1355,'Further physics and validation beyond these implementations',31,True,'#52616E')
proposed=[['Interfacial evaporation','and condensation'],['Compressible waves','and jet formation'],['Species kinetics and','ionization energy'],['Membrane mechanics','and material failure'],['Calibrated optical data','and spectral emission']]
for i,ls in enumerate(proposed):
 x=44+i*389;rect(x,1381,356,101,'#FAFBFC','#8A949D',2,'10 7',9)
 for j,s in enumerate(ls):text(x+178,1421+j*35,s,25,anchor='middle')
text(44,1533,'* Figures 5–6: radial pressure and event-energy updates require correction and convergence tests.',27,True,'#394D5C')
text(44,1575,'Implemented components describe code structure and do not establish physical validation.',27,fill='#52616E')
parts.append('</svg>');svg='\n'.join(parts)
(O/'Figure1_editable.svg').write_text(svg)
# MuPDF converts the vector SVG to a vector PDF and renders a 600 dpi PNG.
s=fitz.open(stream=svg.encode(),filetype='svg');pdf=fitz.open('pdf',s.convert_to_pdf());pdf.save(O/'Figure1_corrected.pdf')
page=pdf[0];page.get_pixmap(dpi=600,alpha=False).save(O/'Figure1_corrected.png')

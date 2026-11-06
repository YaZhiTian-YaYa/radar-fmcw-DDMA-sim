# -*- coding: utf-8 -*-
# 专业重绘 6 张文档图: PIL + 2.5x 超采样 + LANCZOS; 圆角框/菱形/箭头/注记统一设计系统
import os
from PIL import Image, ImageDraw, ImageFont

OUT = r'C:\Users\NJ\AppData\Local\Temp\probe95\docfigs'
SS = 2  # supersample
# ---- 设计系统 ----
C_PROC_F, C_PROC_B = '#E8F1FB', '#4472A8'   # 流程: 淡蓝
C_SE_F,  C_SE_B    = '#FFF3CD', '#B8860B'   # 起止: 淡黄
C_OK_F,  C_OK_B    = '#E6F4EA', '#3D8B4F'   # 输出: 淡绿
C_NO_F,  C_NO_B    = '#FDECEA', '#C0392B'   # 回退: 淡红
C_GATE_F, C_GATE_B = '#FDF2E3', '#C87F0A'   # 判决: 淡橙
C_INK             = '#2B2B2B'
C_GRAY            = '#6E6E6E'
C_ARROW           = '#4A4A4A'
FONT = r'C:\Windows\Fonts\msyh.ttc'
def F(sz): return ImageFont.truetype(FONT, int(sz*SS))

def new_canvas(w, h):
    im = Image.new('RGB', (int(w*SS), int(h*SS)), 'white')
    return im, ImageDraw.Draw(im)

def tlen(d, text, f): return d.textlength(text, font=f)

def wrap(d, text, f, maxw):
    lines, cur = [], ''
    for ch in text:
        if ch == '\n':
            lines.append(cur); cur = ''; continue
        if tlen(d, cur + ch, f) > maxw and cur:
            lines.append(cur); cur = ch
        else:
            cur += ch
    if cur: lines.append(cur)
    return lines

def draw_text_block(d, cx, cy, text, f, maxw, fill=C_INK, lh=None):
    lines = wrap(d, text, f, maxw)
    lh = lh or (f.size * 1.45)
    tot = lh * len(lines)
    y = cy - tot / 2 + lh / 2
    for ln in lines:
        d.text((cx*SS - tlen(d, ln, f)/2, y*SS - f.size*0.72), ln, font=f, fill=fill)
        y += lh
    return tot

def rrect(d, x, y, w, h, r, fill, border, lw=1.6):
    d.rounded_rectangle([x*SS, y*SS, (x+w)*SS, (y+h)*SS], radius=r*SS,
                        fill=fill, outline=border, width=int(lw*SS))

def diamond(d, cx, cy, w, h, fill, border, lw=1.6):
    pts = [(cx*SS, (cy-h/2)*SS), ((cx+w/2)*SS, cy*SS), (cx*SS, (cy+h/2)*SS), ((cx-w/2)*SS, cy*SS)]
    d.polygon(pts, fill=fill, outline=border, width=int(lw*SS))

def box(im_draw, cx, cy, w, h, text, style='proc', shape='round', fsize=13.5, note_side=None):
    d = im_draw
    pal = {'proc': (C_PROC_F, C_PROC_B), 'se': (C_SE_F, C_SE_B), 'ok': (C_OK_F, C_OK_B),
           'no': (C_NO_F, C_NO_B), 'gate': (C_GATE_F, C_GATE_B)}[style]
    if shape == 'diamond':
        diamond(d, cx, cy, w, h, pal[0], pal[1])
    else:
        rrect(d, cx - w/2, cy - h/2, w, h, min(10, h/4), pal[0], pal[1])
    draw_text_block(d, cx, cy, text, F(fsize), w - 22)
    return dict(cx=cx, cy=cy, w=w, h=h)

def arrow(im_draw, p1, p2, label=None, style='solid', lw=2.0, color=C_ARROW):
    d = im_draw
    x1, y1 = p1; x2, y2 = p2
    # 线
    if style == 'dotted':
        dd = ((x2-x1)**2 + (y2-y1)**2) ** .5
        n = max(2, int(dd / 7))
        import math
        for k in range(n):
            t0, t1 = k/n, (k+0.55)/n
            d.line([(x1+(x2-x1)*t0)*SS, (y1+(y2-y1)*t0)*SS,
                    (x1+(x2-x1)*t1)*SS, (y1+(y2-y1)*t1)*SS], fill=color, width=int(lw*SS*0.66))
    else:
        d.line([x1*SS, y1*SS, x2*SS, y2*SS], fill=color, width=int(lw*SS))
    # 箭头头
    import math
    ang = math.atan2(y2-y1, x2-x1)
    L, Wd = 11, 6.5
    p3 = (x2 - L*math.cos(ang) + Wd*math.sin(ang), y2 - L*math.sin(ang) - Wd*math.cos(ang))
    p4 = (x2 - L*math.cos(ang) - Wd*math.sin(ang), y2 - L*math.sin(ang) + Wd*math.cos(ang))
    d.polygon([(x2*SS, y2*SS), (p3[0]*SS, p3[1]*SS), (p4[0]*SS, p4[1]*SS)], fill=color)
    if label:
        f = F(11.5)
        mx, my = (x1+x2)/2, (y1+y2)/2
        d.text((mx*SS - tlen(d, label, f)/2, (my-8)*SS - f.size*0.72), label, font=f, fill=color)

def note(im_draw, boxd, text, dx=26, w=250, align='right'):
    d = im_draw
    x0 = boxd['cx'] + boxd['w']/2 + dx
    y0 = boxd['cy']
    f = F(11)
    lines = wrap(d, text, f, w)
    lh = f.size * 1.4
    y = y0 - lh*len(lines)/2 + lh/2
    for ln in lines:
        d.text((x0*SS, y*SS - f.size*0.72), ln, font=f, fill=C_GRAY)
        y += lh
    arrow(d, (boxd['cx'] + boxd['w']/2 + 3, y0), (x0 - 6, y0), style='dotted', color='#AAAAAA', lw=1.2)

def title(d, cx, y, text, sub=None):
    f = F(21)
    d.text((cx*SS - tlen(d, text, f)/2, y*SS - f.size*0.72), text, font=f, fill=C_INK)
    if sub:
        f2 = F(12.5)
        d.text((cx*SS - tlen(d, sub, f2)/2, (y+24)*SS - f2.size*0.72), sub, font=f2, fill=C_GRAY)

def save(im, name):
    w, h = im.size
    im = im.resize((w//SS, h//SS), Image.LANCZOS)
    im.save(os.path.join(OUT, name))
    print('saved', name)

# ================= 图 1-1 总流程(两遍) =================
im, d = new_canvas(1000, 460)
title(d, 500, 36, 'RCS 标定与估计全流程', '校准遍冻结标定量 → 目标遍按角色逐点估计')
# 校准遍 lane
lane1_y = 84
d.rounded_rectangle([24*SS, lane1_y*SS, 976*SS, (lane1_y+112)*SS], radius=10*SS,
                    fill='#FBFBF7', outline='#D8D4C4', width=int(1.2*SS))
d.text((40*SS, (lane1_y+8)*SS), '校准遍（暗箱，只放角反 3.5 m / −5 dBsm）', font=F(14), fill='#7A6A20')
b1 = box(d, 130, lane1_y+62, 138, 44, '角反测量\n(已知参数)')
b2 = box(d, 296, lane1_y+62, 150, 44, '0°校准表 cCal\n(通道相位掰齐)')
b3 = box(d, 478, lane1_y+62, 158, 44, '噪底+滤形重建\n(档位锚定)')
b4 = box(d, 662, lane1_y+62, 118, 44, 'K0\n(α 域常数)')
b5 = box(d, 810, lane1_y+62, 118, 44, 'K0Coh\n(相参域)', style='ok')
for a, b in [(b1,b2),(b2,b3),(b3,b4),(b4,b5)]:
    arrow(d, (a['cx']+a['w']/2, a['cy']), (b['cx']-b['w']/2, b['cy']))
# 分隔
d.line([(24*SS, 236*SS), (976*SS, 236*SS)], fill='#B9B9B9', width=int(1.4*SS))
for xx in range(24, 976, 26):
    d.line([(xx*SS, 232*SS), ((xx+12)*SS, 240*SS)], fill='#B9B9B9', width=int(1.4*SS))
# 目标遍 lane
lane2_y = 258
d.rounded_rectangle([24*SS, lane2_y*SS, 976*SS, (lane2_y+170)*SS], radius=10*SS,
                    fill='#F7FAF7', outline='#C9D8C9', width=int(1.2*SS))
d.text((40*SS, (lane2_y+8)*SS), '目标遍（冻结标定量，换目标场景）', font=F(14), fill='#2F5F3A')
c1 = box(d, 108, lane2_y+96, 118, 48, '点云点\n五分类角色', style='se')
c2 = box(d, 292, lane2_y+52, 196, 46, '独立点 → 相参图读数\n(峰邻域·插值·俯仰·PLL)', style='proc')
c3 = box(d, 292, lane2_y+142, 196, 46, '配对成员 → α 域联合权重\n(OMP / DML 四类)', style='proc')
c4 = box(d, 536, lane2_y+96, 150, 48, 'σ = SNR + 40lgR\n+ K(支路各配)', style='gate')
c5 = box(d, 760, lane2_y+96, 158, 48, '写回点云 rcsDb\n(+ rcsConf 置信)', style='ok')
arrow(d, (c1['cx']+c1['w']/2, c1['cy']-12), (c2['cx']-c2['w']/2, c2['cy']))
arrow(d, (c1['cx']+c1['w']/2, c1['cy']+14), (c3['cx']-c3['w']/2, c3['cy']))
arrow(d, (c2['cx']+c2['w']/2, c2['cy']), (c4['cx']-c4['w']/2, c4['cy']-10))
arrow(d, (c3['cx']+c3['w']/2, c3['cy']), (c4['cx']-c4['w']/2, c4['cy']+10))
arrow(d, (c4['cx']+c4['w']/2, c4['cy']), (c5['cx']-c5['w']/2, c5['cy']))
arrow(d, (500, lane1_y+112+2), (500, lane2_y-2), label='冻结四量', style='dotted', color='#9A9A9A')
save(im, 'fig_s1_flow.png')

# ================= 图 1-2 五分类路由 =================
im, d = new_canvas(1000, 620)
title(d, 500, 36, '逐点估计的五分类路由', 'splitTag 判角色 → 两条支路各配各的尺')
root = box(d, 175, 110, 190, 56, '点云点\n(角色标记 splitTag)', style='se')
roles = [('独立点 (0)', 'proc'), ('OMP第一目标 (1)', 'proc'), ('OMP第二目标 (2)', 'proc'),
         ('DML第一目标 (3)', 'proc'), ('DML第二目标 (4)', 'proc')]
rb = []
for k, (nm, st) in enumerate(roles):
    y = 78 + k*104
    rb.append(box(d, 430, y, 208, 40, nm, st))
    arrow(d, (root['cx']+root['w']/2, root['cy']), (rb[-1]['cx']-rb[-1]['w']/2, y))
bA = box(d, 742, 176, 216, 74, '相参图读数支路\n（仅独立点）\n尺 = K0Coh', style='ok')
bB = box(d, 742, 390, 216, 74, 'α 域联合权重支路\n（四类配对成员）\n尺 = K0', style='ok')
arrow(d, (rb[0]['cx']+rb[0]['w']/2, rb[0]['cy']), (bA['cx']-bA['w']/2, bA['cy']-16), label='独立点')
arrow(d, (rb[1]['cx']+rb[1]['w']/2, 182), (bB['cx']-bB['w']/2, bB['cy']-18), label='配对成员')
arrow(d, (rb[2]['cx']+rb[2]['w']/2, 286), (bB['cx']-bB['w']/2, bB['cy']-6))
arrow(d, (rb[3]['cx']+rb[3]['w']/2, 390), (bB['cx']-bB['w']/2, bB['cy']+6))
arrow(d, (rb[4]['cx']+rb[4]['w']/2, 494), (bB['cx']-bB['w']/2, bB['cy']+18))
out = box(d, 742, 530, 216, 44, 'σ 写回点云', style='se')
arrow(d, (bA['cx'], bA['cy']+bA['h']/2), (out['cx'], out['cy']-out['h']/2-28))
arrow(d, (bB['cx'], bB['cy']+bB['h']/2), (out['cx'], out['cy']-out['h']/2-2))
note(d, rb[1], 'OMP/DML 的"第一"= 联合解中\n较强者(2026-11-07 修: 原按角度排,\n会出现第二读数高于第一)', dx=20, w=230)
save(im, 'fig_s2_route.png')

# ================= 图 1-3 分离度路由 =================
im, d = new_canvas(1000, 400)
title(d, 500, 36, '同距同速双目标的角度分离度路由')
y0 = 205
d.line([(60*SS, y0*SS), (940*SS, y0*SS)], fill='#3A3A3A', width=int(2.4*SS))
arrow(d, (60, y0), (952, y0))
marks = [(225, '≈2°'), (565, '≈0.7°')]
for xx, lab in marks:
    d.polygon([(xx*SS, (y0-9)*SS), ((xx+9)*SS, y0*SS), (xx*SS, (y0+9)*SS), ((xx-9)*SS, y0*SS)], fill='#3A3A3A')
    f = F(12.5)
    d.text((xx*SS - tlen(d, lab, f)/2, (y0+16)*SS - f.size*0.72), lab, font=f, fill=C_INK)
d.ellipse([(222-5)*SS, (y0-5)*SS, (222+5)*SS, (y0+5)*SS], fill='white', outline='#3A3A3A', width=int(1.6*SS))
d.ellipse([(562-5)*SS, (y0-5)*SS, (562+5)*SS, (y0+5)*SS], fill='white', outline='#3A3A3A', width=int(1.6*SS))
g1 = box(d, 225, 96, 300, 74, 'OMP 区（间隔大）\n残差过 12 dB 门\n→ OMP第一 / OMP第二', style='ok')
g2 = box(d, 565, 96, 264, 74, 'DML 区（0.7°~1.2°）\n双原子联合拟合\n改善>3dB·imprS>1·κ_G≤16', style='gate')
g3 = box(d, 845, 96, 218, 74, '不可分区（<0.7°）\n联合拟合也分不开\n→ 输出融合单点', style='no')
for xx, g in [(225, g1), (565, g2), (845, g3)]:
    arrow(d, (xx, g['cy']+g['h']/2), (xx, y0-6))
n1 = box(d, 225, 300, 300, 44, '第一原子不吸收第二目标', style='proc')
n2 = box(d, 565, 300, 264, 44, '第一原子吸收双方能量，\nOMP 残差塌掉 → 轮到 DML', style='proc')
n3 = box(d, 845, 300, 218, 44, '基底近简并（κ 大）\n单帧信息论极限', style='proc')
for g, n in [(g1, n1), (g2, n2), (g3, n3)]:
    arrow(d, (g['cx'], g['cy']+g['h']/2), (n['cx'], n['cy']-n['h']/2), style='dotted', color='#9A9A9A', lw=1.2)
save(im, 'fig_s3_sep.png')

# ================= 图 1-4 独立点流程 =================
im, d = new_canvas(1000, 1180)
title(d, 500, 34, '独立点 RCS 计算流程', '相参图支路 · 尺 = K0Coh')
steps = [
 ('点云点（角色=独立点）+ 定位 (距离bin, 多普勒bin)', 'se', None),
 ('TXRX 相参图取方位峰：峰 bin ±1 邻域取最大', 'proc', '角栅格 1.4°/bin，目标卡 bin 间有 ±0.7dB 纹波'),
 ('角维三点抛物线定峰（非凹 → 回退原格）', 'proc', '普通读数在 bin 间系统性抖动'),
 ('距离剖面 hannPk3 + 8 根码线多普勒恢复', 'proc', '常规三点抛物线半 bin 处低估 1.18dB（实测 −1.58）'),
 ('俯仰损失补回：θ̂ 取同簇中位，夹 ±15°', 'proc', '单点 θ̂ 有 ±5° 野值（实测多补 +3.5dB）'),
 ('PLL 线损补回 4.343σφ²（芯片参数+自估距离）', 'proc', '角反 3.5m 几乎无损，常数吸不掉（损失 ∝R）'),
 ('SNR = 相参功率读数 − 噪底参考', 'proc', None),
 ('σ = SNR + 40lgR + K0Coh', 'ok', '标定（角反）与估计（目标）同支路同口径，常数相消'),
 ('写回点云 rcsDb（+ rcsConf 置信）', 'se', None)]
y = 108
prev = None
bw, bh, dy = 460, 48, 92
for txt, st, nt in steps:
    b = box(d, 330, y, bw, bh, txt, st)
    if prev: arrow(d, (prev['cx'], prev['cy']+prev['h']/2), (b['cx'], b['cy']-b['h']/2))
    if nt: note(d, b, nt, dx=24, w=250)
    prev = b; y += dy
save(im, 'fig_f_a.png')

# ================= 图 1-5 OMP 流程 =================
im, d = new_canvas(1000, 1280)
title(d, 500, 34, 'OMP 双点 RCS 计算流程', 'α 域联合权重 · 尺 = K0 · 角隔约 ≥2°')
steps = [
 ('检测点（簇）→ 36 通道方位行快拍 x', 'se', None),
 ('原子1 = 导向表插值 @ DBF 方位角', 'proc', None),
 ('LS 扣除第一原子 → 残差匹配曲线', 'proc', None),
 ('@残差峰 − 中位 > 12 dB ？', 'gate', '分离良好时，扣掉强目标后弱目标仍显著'),
 ('第二角度 = 残差峰抛物线亚格定位', 'proc', None),
 ('双原子联合 LS（κ≤4 pinv）', 'proc', '顺序"求一减一"会混入对方泄漏(ρ·s)，\n联合解把互泄漏对消——教科书 OMP 同款'),
 ('按 |W| 排角色：强者为第一（母点），弱者为第二', 'proc', '2026-11-07 修：原母点恒为检测点，\n弱母场景会出现"第二"强于"第一"'),
 ('母点 W1 / 追加点 W2（tag=OMP第一/第二）', 'proc', None),
 ('逐点叶损恢复（距离 hannPk3 + 8 码线均值）', 'proc', '旧口径(插值距离+统计常数)有 ±2dB 幻影，已弃用'),
 ('σ_i = (W_i − 噪底) + 40lgR + K0，两成员各写回', 'ok', None)]
y = 106
prev = None
bw, bh, dy = 470, 50, 98
gate_b = None
for txt, st, nt in steps:
    b = box(d, 330, y, bw, bh, txt, st)
    if prev: arrow(d, (prev['cx'], prev['cy']+prev['h']/2), (b['cx'], b['cy']-b['h']/2))
    if '12 dB' in txt: gate_b = b
    if nt: note(d, b, nt, dx=24, w=252)
    prev = b; y += dy
arrow(d, (gate_b['cx']-gate_b['w']/2, gate_b['cy']), (118, gate_b['cy']), label='否', color='#8A2A20')
nb = box(d, 118, gate_b['cy']+34, 168, 46, '该点按独立点\n走相参图支路', 'no', fsize=12)
save(im, 'fig_f_b.png')

# ================= 图 1-6 DML 流程 =================
im, d = new_canvas(1000, 1340)
title(d, 500, 34, 'DML 双点 RCS 计算流程', 'α 域联合权重 · 尺 = K0 · 角隔 0.7°~1.2°（OMP 检不出的亚波束区）')
steps = [
 ('OMP 未出第二目标的母点 → 36 通道快拍 x', 'se', None),
 ('单原子基线：网格 J1 + 亚bin精化（抛物线角重投影）', 'proc', None),
 ('残差峰候选（主峰 ±8° 内，多起点）', 'proc', None),
 ('交替坐标连续搜索（0.05° 步进 ×3 轮）→ 双角最优', 'proc', '固定一个角精化另一个，轮换到收敛'),
 ('四道门：残差改善>3dB 且 净改善 imprS>1dB\n且 间隔≥0.7° 且 条件数 κ_G≤16', 'gate',
  'imprS 剥掉"亚bin主瓣恢复"被偷算成的改善\n（假 0.3~0.5 / 真等强 3.0）；间隔·条件数门\n拦"两原子塌到强目标一侧"的虚假分离'),
 ('按 |W| 排角色：强者为第一，弱者为第二', 'proc', '2026-11-07 修：原按角度排，等强对分裂噪声\n会使"第二"读得比"第一"高'),
 ('母点方位←强原子角，追加点←弱原子角（tag=DML第一/第二）', 'proc', None),
 ('双原子联合 LS（κ≤4 pinv）→ W1 / W2', 'proc', None),
 ('逐点叶损恢复 → σ_i = (W_i − 噪底) + 40lgR + K0', 'ok', None)]
y = 104
prev = None
bw, bh, dy = 470, 52, 102
gate_b = None
for txt, st, nt in steps:
    b = box(d, 330, y, bw, bh, txt, st)
    if prev: arrow(d, (prev['cx'], prev['cy']+prev['h']/2), (b['cx'], b['cy']-b['h']/2))
    if '四道门' in txt: gate_b = b
    if nt: note(d, b, nt, dx=24, w=258)
    prev = b; y += dy
arrow(d, (gate_b['cx']-gate_b['w']/2, gate_b['cy']), (118, gate_b['cy']), label='否', color='#8A2A20')
nb = box(d, 116, gate_b['cy']+44, 172, 60, '按设计输出融合单点\n（独立支路读数\n≈ 两目标相干和）', 'no', fsize=11.5)
arrow(d, (nb['cx'], nb['cy']+nb['h']/2), (nb['cx'], nb['cy']+96), style='dotted', color='#9A9A9A')
d.text(((116-84)*SS, (nb['cy']+112)*SS), '单帧口径终点行为（车载高速不用多帧）', font=F(11), fill=C_GRAY)
save(im, 'fig_f_c.png')
print('ALL FIGS DONE')

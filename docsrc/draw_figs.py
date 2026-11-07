# -*- coding: utf-8 -*-
# 专业图示引擎 v2: 文字测量→框自适应; 正交总线布线(避让所有框); 箭头止于框边
import os, math
from PIL import Image, ImageDraw, ImageFont

OUT = r'C:\Users\NJ\AppData\Local\Temp\probe95\docfigs'
SS = 2
C_PROC_F, C_PROC_B = '#E8F1FB', '#4472A8'
C_SE_F,  C_SE_B    = '#FFF3CD', '#B8860B'
C_OK_F,  C_OK_B    = '#E6F4EA', '#3D8B4F'
C_NO_F,  C_NO_B    = '#FDECEA', '#C0392B'
C_GATE_F, C_GATE_B = '#FDF2E3', '#C87F0A'
C_INK, C_GRAY, C_ARROW = '#2B2B2B', '#6E6E6E', '#4A4A4A'
PAL = {'proc': (C_PROC_F, C_PROC_B), 'se': (C_SE_F, C_SE_B), 'ok': (C_OK_F, C_OK_B),
       'no': (C_NO_F, C_NO_B), 'gate': (C_GATE_F, C_GATE_B)}
FONT = r'C:\Windows\Fonts\msyh.ttc'
def F(sz): return ImageFont.truetype(FONT, int(sz*SS))

class Fig:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.im = Image.new('RGB', (w*SS, h*SS), 'white')
        self.d = ImageDraw.Draw(self.im)
        self.boxes = []

    # ---- 文字测量 ----
    def wrap(self, text, f, maxw):
        lines, cur = [], ''
        for ch in text:
            if ch == '\n':
                lines.append(cur); cur = ''; continue
            w = self.d.textlength(cur+ch, font=f)
            if w > maxw and cur: lines.append(cur); cur = ch
            else: cur += ch
        if cur or not lines: lines.append(cur)
        return lines

    def text_h(self, lines, f): return len(lines) * f.size / SS * 1.5

    # ---- 框: 自适应 ----
    def box(self, cx, cy, text, style='proc', w=None, min_h=40, fsize=13.5, minh_only=True):
        f = F(fsize)
        # 自适应宽: 先按期望宽换行, 若只有 1 行且很短 → 收窄; 行数过多 → 加宽重排
        w0 = w or 420
        for _ in range(6):
            lines = self.wrap(text, f, (w0-26)*SS)
            if len(lines) <= 3 or w is not None and w0 >= w: break
            w0 += 40
        maxlw = max(self.d.textlength(ln, font=f) for ln in lines) / SS
        bw = w0 if w else min(w0, maxlw + 30)
        bh = max(min_h, self.text_h(lines, f) + 22)
        fill, border = PAL[style]
        x0, y0 = cx - bw/2, cy - bh/2
        self.d.rounded_rectangle([x0*SS, y0*SS, (x0+bw)*SS, (y0+bh)*SS],
                                 radius=min(10, bh/4)*SS, fill=fill, outline=border, width=int(1.6*SS))
        # 居中画行
        lh = f.size/SS*1.5
        ty = cy - lh*len(lines)/2 + lh/2
        for ln in lines:
            self.d.text((cx*SS - self.d.textlength(ln, font=f)/2, ty*SS - f.size*0.74), ln, font=f, fill=C_INK)
            ty += lh
        b = dict(cx=cx, cy=cy, w=bw, h=bh, x0=x0, y0=y0, x1=x0+bw, y1=y0+bh)
        self.boxes.append(b)
        return b

    # ---- 边缘锚点 & 无碰撞检查 ----
    def edge_pt(self, b, side):
        return {'l': (b['x0'], b['cy']), 'r': (b['x1'], b['cy']),
                't': (b['cx'], b['y0']), 'b': (b['cx'], b['y1'])}[side]

    def seg_clear(self, p1, p2, skip=()):
        for b in self.boxes:
            if b in skip: continue
            if self._seg_rect(p1, p2, b): return False
        return True

    def _seg_rect(self, p1, p2, b):
        x1, y1 = p1; x2, y2 = p2
        rx0, ry0, rx1, ry1 = b['x0'], b['y0'], b['x1'], b['y1']
        if x1 == x2:  # 竖直: 严格内部(x 不等于边线)
            return rx0 < x1 < rx1 and max(y1,y2) > ry0 and min(y1,y2) < ry1
        if y1 == y2:  # 水平: 严格内部
            return ry0 < y1 < ry1 and max(x1,x2) > rx0 and min(x1,x2) < rx1
        return True   # 斜线按保守处理(本引擎不产生斜线)

    # ---- 直线箭头(带碰撞检查) ----
    def arrow(self, p1, p2, label=None, color=C_ARROW, lw=2.0, check=True, skip=()):
        if check and not self.seg_clear(p1, p2, skip):
            raise RuntimeError('直线箭头穿框: %s -> %s' % (p1, p2))
        self._line(p1, p2, color, lw, solid=True)
        self._head(p1, p2, color)
        if label:
            f = F(11.5)
            mx, my = (p1[0]+p2[0])/2, (p1[1]+p2[1])/2
            self.d.text((mx*SS - self.d.textlength(label, font=f)/2, (my-9)*SS - f.size*0.74), label, font=f, fill=color)

    # ---- 正交布线: 出边 → 总线 → 入边 (Z 形) ----
    def ortho(self, a, b, aside, bside, bus=None, label=None, color=C_ARROW, lw=2.0):
        pa = self.edge_pt(a, aside); pb = self.edge_pt(b, bside)
        if bus is None:
            bus = (pa[0] + pb[0])/2 if aside in 'lr' else (pa[1] + pb[1])/2
        if aside in 'lr':
            segs = [pa, (bus, pa[1]), (bus, pb[1]), pb]
        else:
            segs = [pa, (pa[0], bus), (pb[0], bus), pb]
        for s0, s1 in zip(segs[:-1], segs[1:]):
            if not self.seg_clear(s0, s1, skip=(a, b)):
                raise RuntimeError('正交连线穿框: %s -> %s' % (s0, s1))
        for s, e in zip(segs[:-1], segs[1:-1]):
            self._line(s, e, color, lw, solid=True)
        self._line(segs[-2], segs[-1], color, lw, solid=True)
        self._head(segs[-2], pb, color)
        if label:
            f = F(11.5)
            lx = (pa[0]+bus)/2 if aside in 'lr' else pa[0]
            ly = pa[1] - 14 if aside in 'lr' else (pa[1]+bus)/2
            self.d.text((lx*SS - self.d.textlength(label, font=f)/2, ly*SS - f.size*0.74), label, font=f, fill=color)

    def _line(self, p1, p2, color, lw, solid=True):
        if solid:
            self.d.line([p1[0]*SS, p1[1]*SS, p2[0]*SS, p2[1]*SS], fill=color, width=int(lw*SS))
        else:
            dd = math.hypot(p2[0]-p1[0], p2[1]-p1[1]); n = max(2, int(dd/7))
            for k in range(n):
                t0, t1 = k/n, (k+0.55)/n
                self.d.line([(p1[0]+(p2[0]-p1[0])*t0)*SS, (p1[1]+(p2[1]-p1[1])*t0)*SS,
                             (p1[0]+(p2[0]-p1[0])*t1)*SS, (p1[1]+(p2[1]-p1[1])*t1)*SS],
                            fill=color, width=int(lw*SS*0.66))

    def _head(self, p1, p2, color):
        ang = math.atan2(p2[1]-p1[1], p2[0]-p1[0])
        L, W = 11, 6.5
        p3 = (p2[0]-L*math.cos(ang)+W*math.sin(ang), p2[1]-L*math.sin(ang)-W*math.cos(ang))
        p4 = (p2[0]-L*math.cos(ang)-W*math.sin(ang), p2[1]-L*math.sin(ang)+W*math.cos(ang))
        self.d.polygon([(p2[0]*SS, p2[1]*SS), (p3[0]*SS, p3[1]*SS), (p4[0]*SS, p4[1]*SS)], fill=color)

    # ---- 注记(右侧, 虚线连到框边) ----
    def note(self, b, text, dx=22, w=252):
        x0 = b['x1'] + dx
        f = F(11)
        lines = self.wrap(text, f, w*SS)
        lh = f.size/SS*1.42
        y = b['cy'] - lh*len(lines)/2 + lh/2
        for ln in lines:
            self.d.text((x0*SS, y*SS - f.size*0.74), ln, font=f, fill=C_GRAY)
            y += lh
        self._line((b['x1']+3, b['cy']), (x0-6, b['cy']), '#AAAAAA', 1.2, solid=False)

    def title(self, text, sub=None):
        f = F(21)
        self.d.text((self.w/2*SS - self.d.textlength(text, font=f)/2, 36*SS - f.size*0.74), text, font=f, fill=C_INK)
        if sub:
            f2 = F(12.5)
            self.d.text((self.w/2*SS - self.d.textlength(sub, font=f2)/2, 62*SS - f2.size*0.74), sub, font=f2, fill=C_GRAY)

    def save(self, name):
        im = self.im.resize((self.w*SS//SS, self.h*SS//SS), Image.LANCZOS)
        im.save(os.path.join(OUT, name))
        print('saved', name)

# ================= 图 1-1 =================
g = Fig(1000, 470)
g.title('RCS 标定与估计全流程', '校准遍冻结标定量 → 目标遍按角色逐点估计')
g.d.rounded_rectangle([24*SS, 84*SS, 976*SS, 196*SS], radius=10*SS, fill='#FBFBF7', outline='#D8D4C4', width=int(1.2*SS))
g.d.text((40*SS, 92*SS), '校准遍（暗箱，只放角反 3.5 m / −5 dBsm）', font=F(14), fill='#7A6A20')
b1 = g.box(130, 146, '角反测量\n(已知参数)', 'se', w=150)
b2 = g.box(296, 146, '0°校准表 cCal\n(通道相位校准对齐)', 'proc', w=168)
b3 = g.box(482, 146, '噪底+滤形重建\n(档位锚定)', 'proc', w=176)
b4 = g.box(668, 146, 'K0\n(α 域常数)', 'proc', w=118)
b5 = g.box(830, 146, 'K0Coh\n(相参域)', 'ok', w=112)
for a, b in [(b1,b2),(b2,b3),(b3,b4),(b4,b5)]:
    g.arrow((a['x1'], a['cy']), (b['x0'], b['cy']))
for xx in range(24, 976, 26):
    g.d.line([(xx*SS, 224*SS), ((xx+12)*SS, 232*SS)], fill='#B9B9B9', width=int(1.4*SS))
g.d.line([(24*SS, 228*SS), (976*SS, 228*SS)], fill='#B9B9B9', width=int(1.4*SS))
g.d.rounded_rectangle([24*SS, 250*SS, 976*SS, 440*SS], radius=10*SS, fill='#F7FAF7', outline='#C9D8C9', width=int(1.2*SS))
g.d.text((40*SS, 258*SS), '目标遍（冻结标定量，换目标场景）', font=F(14), fill='#2F5F3A')
c1 = g.box(112, 340, '点云点\n五分类角色', 'se', w=132)
c2 = g.box(300, 296, '独立点 → 相参图读数\n峰邻域·插值·俯仰·PLL', 'proc', w=216)
c3 = g.box(300, 386, '配对成员 → α 域联合权重\n(OMP / DML 四类)', 'proc', w=216)
c4 = g.box(560, 340, 'σ = SNR + 40lgR\n+ K（按支路标定）', 'gate', w=168)
c5 = g.box(790, 340, '写回点云 rcsDb\n(+ rcsConf 置信)', 'ok', w=168)
g.ortho(c1, c2, 'r', 'l', bus=224)
g.ortho(c1, c3, 'r', 'l', bus=224)
g.ortho(c2, c4, 'r', 'l', bus=470)
g.ortho(c3, c4, 'r', 'l', bus=470)
g.arrow((c4['x1'], c4['cy']), (c5['x0'], c5['cy']))
g._line((500, 198), (500, 248), '#9A9A9A', 1.4, solid=False)
g._head((500, 198), (500, 248), '#9A9A9A')
fz = F(11.5)
g.d.text((512*SS, 216*SS - fz.size*0.74), '冻结四量', font=fz, fill='#9A9A9A')
g.save('fig_s1_flow.png')

# ================= 图 1-2 =================
g = Fig(1000, 640)
g.title('逐点估计的五分类路由', 'splitTag 判角色 → 两条支路分别标定')
root = g.box(160, 116, '点云点\n(角色标记 splitTag)', 'se', w=192)
roles = ['独立点 (0)', 'OMP第一目标 (1)', 'OMP第二目标 (2)', 'DML第一目标 (3)', 'DML第二目标 (4)']
rb = []
for k, nm in enumerate(roles):
    rb.append(g.box(430, 76 + k*106, nm, 'proc', w=210))
    g.ortho(root, rb[-1], 'r', 'l', bus=316)
bA = g.box(756, 176, '相参图读数支路\n（仅独立点）\n尺 = K0Coh', 'ok', w=210)
bB = g.box(756, 400, 'α 域联合权重支路\n（四类配对成员）\n尺 = K0', 'ok', w=210)
BUS2 = 596
g.ortho(rb[0], bA, 'r', 'l', bus=BUS2)
for k in [1, 2, 3, 4]:
    g.ortho(rb[k], bB, 'r', 'l', bus=BUS2)
outb = g.box(716, 560, 'σ 写回点云', 'se', w=170)
segsB = [((bB['cx'], bB['y1']), (bB['cx'], (bB['y1']+outb['y0'])/2)),
         ((bB['cx'], (bB['y1']+outb['y0'])/2), (outb['cx'], (bB['y1']+outb['y0'])/2)),
         ((outb['cx'], (bB['y1']+outb['y0'])/2), (outb['cx'], outb['y0']))]
for s0, s1 in segsB:
    if not g.seg_clear(s0, s1, skip=(bB, outb)):
        raise RuntimeError('bB 连线穿框: %s' % (s0,))
for s0, s1 in segsB[:-1]:
    g._line(s0, s1, C_ARROW, 2.0)
g._line(segsB[-1][0], segsB[-1][1], C_ARROW, 2.0)
g._head(segsB[-1][0], segsB[-1][1], C_ARROW)
CH = 905
segsA = [((bA['cx'], bA['y1']), (bA['cx'], bA['y1']+24)),
         ((bA['cx'], bA['y1']+24), (CH, bA['y1']+24)),
         ((CH, bA['y1']+24), (CH, outb['cy'])),
         ((CH, outb['cy']), (outb['x1']+8, outb['cy']))]
for s0, s1 in segsA:
    if not g.seg_clear(s0, s1, skip=(bA, outb)):
        raise RuntimeError('绕行通道穿框: %s' % (s0,))
for s0, s1 in segsA[:-1]:
    g._line(s0, s1, C_ARROW, 2.0)
g._line(segsA[-1][0], segsA[-1][1], C_ARROW, 2.0)
g._head(segsA[-1][0], segsA[-1][1], C_ARROW)
# 角色排序注记(连到 OMP第一 那行)
g.save('fig_s2_route.png')

# ================= 图 1-3 =================
g = Fig(1000, 420)
g.title('同距同速双目标的角度分离度路由')
y0 = 212
g._line((60, y0), (952, y0), '#3A3A3A', 2.4)
g._head((60, y0), (952, y0), '#3A3A3A')
for xx, lab in [(225, '≈2°'), (565, '≈0.7°')]:
    g.d.polygon([(xx*SS, (y0-9)*SS), ((xx+9)*SS, y0*SS), (xx*SS, (y0+9)*SS), ((xx-9)*SS, y0*SS)], fill='#3A3A3A')
    f = F(12.5)
    g.d.text((xx*SS - g.d.textlength(lab, font=f)/2, (y0+16)*SS - f.size*0.74), lab, font=f, fill=C_INK)
    g.d.ellipse([(xx-5)*SS, (y0-5)*SS, (xx+5)*SS, (y0+5)*SS], fill='white', outline='#3A3A3A', width=int(1.6*SS))
g1 = g.box(225, 100, 'OMP 区（间隔大）\n残差过 12 dB 门\n→ OMP第一 / OMP第二', 'ok', w=300)
g2 = g.box(565, 100, 'DML 区（0.7°~1.2°）\n双原子联合拟合\n改善>3dB·imprS>1·κ_G≤16', 'gate', w=270)
g3 = g.box(848, 100, '不可分区（<0.7°）\n联合拟合也分不开\n→ 输出融合单点', 'no', w=216)
for xx, gb in [(225, g1), (565, g2), (848, g3)]:
    g.arrow((xx, gb['y1']), (xx, y0-6), check=False)
n1 = g.box(225, 316, '第一原子不吸收第二目标', 'proc', w=280)
n2 = g.box(565, 316, '第一原子吸收双方能量，OMP 残差显著降低', 'proc', w=280)
n3 = g.box(848, 316, '基底近简并（κ 大）单帧信息论极限', 'proc', w=280)
for gb, nb in [(g1, n1), (g2, n2), (g3, n3)]:
    g._line((gb['cx'], gb['y1']+2), (nb['cx'], nb['y0']-2), '#9A9A9A', 1.2, solid=False)
g.save('fig_s3_sep.png')

# ================= 图 1-4 独立点 =================
g = Fig(1000, 1210)
g.title('独立点 RCS 计算流程', '相参图支路 · 尺 = K0Coh')
steps = [
 ('点云点（角色=独立点）+ 定位 (距离bin, 多普勒bin)', 'se', None),
 ('TXRX 相参图取方位峰：峰 bin ±1 邻域取最大', 'proc', '角栅格 1.4°/bin，目标位于单元之间\n存在 ±0.7dB 纹波'),
 ('角维三点抛物线定峰（非凹 → 回退原格）', 'proc', '普通读数在 bin 间系统性抖动'),
 ('距离剖面 hannPk3 + 8 根码线多普勒恢复', 'proc', '常规三点抛物线半单元偏移处低估\n1.18dB（实测 −1.58）'),
 ('俯仰损失补回：θ̂ 取同簇中位，夹 ±15°', 'proc', '单点 θ̂ 存在 ±5° 级野值\n（实测出现过补偿 +3.5dB）'),
 ('PLL 线损补回 4.343σφ²（芯片参数+自估距离）', 'proc', '角反 3.5m 近乎无损，常数无法吸收\n（损失 ∝ 距离）'),
 ('SNR = 相参功率读数 − 噪底参考', 'proc', None),
 ('σ = SNR + 40lgR + K0Coh', 'ok', '标定（角反）与估计（目标）同支路\n同口径，标定误差在比值中相消'),
 ('写回点云 rcsDb（+ rcsConf 置信）', 'se', None)]
y, prev, gap = 106, None, 22
BW = 470
for txt, st, nt in steps:
    b = g.box(320, y, txt, st, w=BW)
    if prev: g.arrow((prev['cx'], prev['y1']), (b['cx'], b['y0']))
    if nt: g.note(b, nt, dx=22, w=248)
    prev = b; y = b['y1'] + gap + (b['h'] if False else 0) + 34
g.save('fig_f_a.png')

# ================= 图 1-5 OMP =================
g = Fig(1020, 1310)
g.title('OMP 双点 RCS 计算流程', 'α 域联合权重 · 尺 = K0 · 角隔约 ≥2°')
steps = [
 ('检测点（簇）→ 36 通道方位行快拍 x', 'se', None),
 ('原子1 = 导向表插值 @ DBF 方位角', 'proc', None),
 ('LS 扣除第一原子 → 残差匹配曲线', 'proc', None),
 ('残差峰 − 中位 > 12 dB ？', 'gate', '分离良好时，扣掉强目标后\n弱目标仍显著'),
 ('第二角度 = 残差峰抛物线亚格定位', 'proc', None),
 ('双原子联合 LS（κ≤4 pinv）', 'proc', '顺序估计会混入对方泄漏\n(ρ·s)，联合求解将互泄漏对消\n——与标准 OMP 实现一致'),
 ('按 |W| 排角色：强者为第一（母点）', 'proc', '2026-11-07 修：原实现母点恒为检测\n点，弱目标先被检测时第二目标读数会高于第一目标'),
 ('母点 W1 / 追加点 W2（tag=OMP第一/第二）', 'proc', None),
 ('逐点叶损恢复（距离 hannPk3 + 8 码线均值）', 'proc', '旧方法（插值距离+统计常数）存在\n±2dB 误差，已弃用'),
 ('σ_i = (W_i − 噪底) + 40lgR + K0，两成员各写回', 'ok', None)]
y, prev, gate_b = 104, None, None
BW = 452
for txt, st, nt in steps:
    b = g.box(462, y, txt, st, w=BW)
    if prev: g.arrow((prev['cx'], prev['y1']), (b['cx'], b['y0']))
    if '12 dB' in txt: gate_b = b
    if nt: g.note(b, nt, dx=22, w=236)
    prev = b; y = b['y1'] + 58
nb = g.box(122, gate_b['cy'], '该点按独立点\n走相参谱支路', 'no', w=158, fsize=12)
g.arrow((gate_b['x0'], gate_b['cy']), (nb['x1'], nb['cy']), label='否', skip=(gate_b, nb), color='#8A2A20')
g.save('fig_f_b.png')

# ================= 图 1-6 DML =================
g = Fig(1020, 1400)
g.title('DML 双点 RCS 计算流程', 'α 域联合权重 · 尺 = K0 · 角隔 0.7°~1.2°（OMP 检不出的亚波束区）')
steps = [
 ('OMP 未出第二目标的母点 → 36 通道快拍 x', 'se', None),
 ('单原子基线：网格 J1 + 亚bin精化（抛物线角重投影）', 'proc', None),
 ('残差峰候选（主峰 ±8° 内，多起点）', 'proc', None),
 ('交替坐标连续搜索（0.05° 步进 ×3 轮）→ 双角最优', 'proc', '固定一角精化另一角，\n交替迭代至收敛'),
 ('四道门：残差改善>3dB 且 净改善 imprS>1dB\n且 间隔≥0.7° 且 条件数 κ_G≤16', 'gate',
  'imprS 剥掉"亚bin主瓣恢复"被错误计入\n的改善（假 0.3~0.5 / 真等强 3.0）；\n间隔·条件数门抑制双原子坍缩至强目标\n一侧的虚假分离'),
 ('按 |W| 排角色：强者为第一，弱者为第二', 'proc', '2026-11-07 修：原按角度序分配，等强\n对的分裂噪声曾导致第二目标读数高于第一目标'),
 ('母点方位←强原子角，追加点←弱原子角', 'proc', None),
 ('双原子联合 LS（κ≤4 pinv）→ W1 / W2', 'proc', None),
 ('逐点叶损恢复 → σ_i = (W_i − 噪底) + 40lgR + K0', 'ok', None)]
y, prev, gate_b = 102, None, None
BW = 452
for txt, st, nt in steps:
    b = g.box(462, y, txt, st, w=BW)
    if prev: g.arrow((prev['cx'], prev['y1']), (b['cx'], b['y0']))
    if '四道门' in txt: gate_b = b
    if nt: g.note(b, nt, dx=22, w=232)
    prev = b; y = b['y1'] + 60
nb = g.box(126, gate_b['cy'], '按设计输出融合单点\n（独立支路读数\n≈ 两目标相干和）', 'no', w=166, fsize=11.5)
g.arrow((gate_b['x0'], gate_b['cy']), (nb['x1'], nb['cy']), label='否', skip=(gate_b, nb), color='#8A2A20')
g._line((nb['cx'], nb['y1']), (nb['cx'], nb['y1']+42), '#9A9A9A', 1.2, solid=False)
fn = F(11)
g.d.text((nb['cx']*SS - g.d.textlength('单帧体制下的最终行为（车载高速场景不采用多帧）', font=fn)/2,
          (nb['y1']+48)*SS - fn.size*0.74), '单帧体制下的最终行为（车载高速场景不采用多帧）', font=fn, fill=C_GRAY)
g.save('fig_f_c.png')
print('ALL FIGS DONE v2')

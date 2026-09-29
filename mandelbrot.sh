#!/usr/bin/env bash
# mandelbrot.sh - full-screen animated ASCII Mandelbrot zoom
#
# Usage: ./mandelbrot.sh [-x fractal] [-m] [-t target] [-z secs] [-a] [-c|-o] [-l steps] [-s n] [-r] [-w] [-f size]
#   -x N        fractal: 1 = Mandelbrot (default), 2 = Multibrot z^3,
#               3 = Multibrot z^4, 4 = Burning Ship, 5 = Tricorn,
#               6 = Julia (morphing), 0 = all of them, one after another
#   -m          monochrome (no ANSI colors)
#   -t N        start at the fractal's Nth zoom target (default 1). Mandelbrot:
#               1 Seahorse spiral, 2 Elephant Valley, 3 Double spiral,
#               4 Quad spiral arms, 5 Antenna feather, 6 Dendrite branch;
#               the other fractals have 3 each. After each zoom the next
#               target follows. The status bar shows which one is playing.
#   -z SECS     zoom speed: seconds per 10x zoom (default 6, higher = slower)
#   -a          autopilot: wander around, steering towards detail and away
#               from the black inside (default: steady straight zoom)
#   -c          cycle (default): zoom into each pattern for 1000 steps (frames),
#               rewind back out, then move on to the next pattern, looping all 6
#   -o          one long zoom per pattern, down to the precision limit
#   -l STEPS    steps per zoom-in when cycling (default 1000)
#   -s N        smoothing: N x N points per character, 1-4 (default 3;
#               1 = fastest, 3-4 = smoothest but slower)
#   -r          use the long 57-character ramp instead of the default
#               15-character " .,:-~=+*o#%&8@" (finer but busier shading)
#   -w          open a new fullscreen foot window with a small font
#               (default: render in the current terminal, at its size)
#   -f SIZE     font size for -w (default 5)
# Quit with Ctrl+C or q.

COLOR=1
TARGET=1
FRACTAL=1
FONTSIZE=5
NEWWIN=0
ZSECS=6
AUTOPILOT=0
SS=3
RAMP=" .,:-~=+*o#%&8@"
STEPS=1000
ARGS=("$@")
while getopts "x:mt:z:acol:s:rf:wnh" opt; do
  case $opt in
    x) FRACTAL=$OPTARG ;;
    m) COLOR=0 ;;
    t) TARGET=$OPTARG ;;
    z) ZSECS=$OPTARG ;;
    a) AUTOPILOT=1 ;;
    c) (( STEPS )) || STEPS=1000 ;;
    o) STEPS=0 ;;
    l) STEPS=$OPTARG ;;
    s) SS=$OPTARG ;;
    r) RAMP="" ;;
    f) FONTSIZE=$OPTARG ;;
    w) NEWWIN=1 ;;
    n) NEWWIN=0 ;;
    *) sed -n '2,28p' "$0"; exit 0 ;;
  esac
done

# ---------------------------------------------------------------------------
# Fast renderer: a small C program, compiled once and cached. It renders on
# all CPU cores (OpenMP). If no C compiler is found, the awk renderer is used.
# ---------------------------------------------------------------------------
C_SRC=$(cat <<'EOF'
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* Character ramp ordered from sparse (far from set) to dense (at the edge) */
static const char default_ramp[] =
    " .'`^\",:;l!i><~+_-?][}{1)(|\\/tfjrxnuvcz0mwqpdbkhao*#&8%@$";

#define NQ 256   /* resolution of the per-frame escape-time histogram */

static float qtab[NQ + 1]; /* escape-time quantiles of the previous frame */
static int have_q = 0, qready = 0;

/* Fractal types */
enum { F_MANDEL = 1, F_MULTI3, F_MULTI4, F_SHIP, F_TRICORN, F_JULIA };
static int ftype = F_MANDEL;
static double lnd = 0.6931471805599453; /* log(degree) for smooth colouring */
static double jr = 0, ji = 0;            /* Julia parameter */

/* Iterate the fractal's formula. Returns 1 and the smooth iteration count
   if the point escapes, 0 if it is (probably) inside the set. */
static inline int escapes(double cr, double ci, int maxit, double *nu_out) {
    double zr = 0, zi = 0;
    if (ftype == F_JULIA) {
        zr = cr; zi = ci; cr = jr; ci = ji;
    } else if (ftype == F_MANDEL) {
        /* quick reject: main cardioid and period-2 bulb */
        double q = (cr - 0.25) * (cr - 0.25) + ci * ci;
        if (q * (q + cr - 0.25) < 0.25 * ci * ci || (cr + 1) * (cr + 1) + ci * ci < 0.0625)
            return 0;
    }
    double m2 = zr * zr + zi * zi, pr = zr, pi = zi, t;
    int it = 0, check = 8;
    while (m2 < 256) {
        if (it >= maxit) return 0;
        switch (ftype) {
        case F_MULTI3: { /* z^3 + c */
            double r2 = zr * zr, i2 = zi * zi;
            t = zr * (r2 - 3 * i2) + cr;
            zi = zi * (3 * r2 - i2) + ci;
            zr = t;
        } break;
        case F_MULTI4: { /* z^4 + c */
            double sr = zr * zr - zi * zi, si = 2 * zr * zi;
            t = sr * sr - si * si + cr;
            zi = 2 * sr * si + ci;
            zr = t;
        } break;
        case F_SHIP: { /* (|x| + i|y|)^2 + c */
            double ar = fabs(zr), ai = fabs(zi);
            t = ar * ar - ai * ai + cr;
            zi = 2 * ar * ai + ci;
            zr = t;
        } break;
        case F_TRICORN: /* conj(z)^2 + c */
            t = zr * zr - zi * zi + cr;
            zi = -2 * zr * zi + ci;
            zr = t;
            break;
        default: /* Mandelbrot and Julia: z^2 + c */
            t = zr * zr - zi * zi + cr;
            zi = 2 * zr * zi + ci;
            zr = t;
        }
        m2 = zr * zr + zi * zi;
        it++;
        /* periodicity check: orbit returning to a saved point means it is
           caught in a cycle, i.e. inside the set */
        if (fabs(zr - pr) < 1e-16 && fabs(zi - pi) < 1e-16) return 0;
        if (it == check) { pr = zr; pi = zi; check *= 2; }
    }
    double nu = it + 1 - log(log(m2) / 2 / 0.6931471805599453) / lnd;
    *nu_out = nu < 1 ? 1 : nu;
    return 1;
}

/* Position 0..1 on the character ramp for a smooth iteration count */
static inline double ramp_pos(double nu, double logmax) {
    double t;
    if (have_q) {
        /* histogram equalisation: rank of nu among last frame's values,
           so every frame uses the whole ramp */
        int lo = 0, hi = NQ;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (qtab[mid] <= nu) lo = mid; else hi = mid;
        }
        double span = qtab[hi] - qtab[lo];
        t = (lo + (span > 0 ? (nu - qtab[lo]) / span : 0)) / NQ;
    } else {
        t = log(nu) / logmax;
    }
    return t < 0 ? 0 : t > 1 ? 1 : t;
}

static int cmpf(const void *a, const void *b) {
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv) {
    if (argc < 6) return 1;
    int W = atoi(argv[1]), H = atoi(argv[2]) - 1; /* last line = status bar */
    double CX = atof(argv[3]), CY = atof(argv[4]);
    int COLOR = atoi(argv[5]);
    double zsecs = argc > 6 ? atof(argv[6]) : 5; /* seconds per 10x zoom */
    int autopilot = argc > 7 ? atoi(argv[7]) : 1;
    int SS = argc > 8 ? atoi(argv[8]) : 2; /* SS x SS samples per character */
    if (SS < 1) SS = 1;
    if (SS > 4) SS = 4;
    if (zsecs < 0.1) zsecs = 0.1;
    const char *ramp = argc > 9 && argv[9][0] ? argv[9] : default_ramp;
    int steps = argc > 10 ? atoi(argv[10]) : 0; /* >0: zoom N steps, rewind, exit */
    ftype = argc > 11 ? atoi(argv[11]) : F_MANDEL;
    if (ftype < F_MANDEL || ftype > F_JULIA) ftype = F_MANDEL;
    lnd = log(ftype == F_MULTI3 ? 3.0 : ftype == F_MULTI4 ? 4.0 : 2.0);
    const char *label = argc > 12 ? argv[12] : "";
    int n = (int)strlen(ramp);
    const double aspect = 2.1; /* terminal cells are ~2x taller than wide */
    const double frame_ns = 1e9 / 30; /* cap at 30 fps */
    if (W < 2 || H < 1) return 1;

    /* Smooth 256-colour palette (cube colours) */
    char pal[64][16];
    int pal_len[64];
    for (int i = 0; i < 64; i++) {
        double t = i / 64.0;
        int r = (int)(5 * (0.5 + 0.5 * sin(6.2832 * (t + 0.00))));
        int g = (int)(5 * (0.5 + 0.5 * sin(6.2832 * (t + 0.33))));
        int b = (int)(5 * (0.5 + 0.5 * sin(6.2832 * (t + 0.67))));
        pal_len[i] = snprintf(pal[i], 16, "\033[38;5;%dm", 16 + 36 * r + 6 * g + b);
    }
    static const char inside_col[] = "\033[38;5;232m";

    size_t rowcap = (size_t)W * 13 + 2;
    char *rows = malloc(rowcap * H);
    int *rowlen = malloc(sizeof(int) * H);
    char *out = malloc(rowcap * H + 1024);
    float *nubuf = malloc(sizeof(float) * W * H); /* smooth iter, -1 = inside */
    float *tbuf = malloc(sizeof(float) * W * H);  /* ramp position 0..1 */
    float *samp = malloc(sizeof(float) * 8192);
    if (!rows || !rowlen || !out || !nubuf || !tbuf || !samp) return 1;

    const double scale0 = (ftype == F_TRICORN ? 4.0 : ftype == F_MULTI3 || ftype == F_MULTI4 || ftype == F_JULIA ? 3.2 : 3.5) / W;
    double scale = scale0, zoom = 1;
    int frame = 0;
    struct timespec t0, t1;
    double fps = 0, dt = 1.0 / 30;
    double TR = CX, TI = CY;     /* autopilot target (complex plane) */
    double black = 0, zspeed = 1, zsmooth = 0, elapsed = 0;

    /* cycle mode state: L = zoom depth in decades, v = its velocity */
    int step = 0, rewinding = 0;
    double L = 0, v = 0;
    const double dz = 1.0 / (30 * zsecs);     /* decades per step */
    /* max zoom depth: double precision limit; a morphing Julia set has no
       fixed point to dive into, so it only zooms in a little */
    const double Lmax = ftype == F_JULIA ? 0.15 : log10(scale0 / 5e-15);

    while (steps || scale > 2e-15) { /* stop at double precision limit */
        clock_gettime(CLOCK_MONOTONIC, &t0);
        if (ftype == F_JULIA) {
            /* slowly move c along a path just outside the Mandelbrot set's
               main cardioid, so the Julia set keeps morphing */
            /* driven by the wall clock, so the morph carries on across cycles */
            struct timespec now;
            clock_gettime(CLOCK_REALTIME, &now);
            double a = fmod(now.tv_sec % 100000 + now.tv_nsec / 1e9, 1e5) * 0.08, wr = cos(a), wi = sin(a);
            jr = 1.02 * (wr / 2 - (wr * wr - wi * wi) / 4);
            ji = 1.02 * (wi / 2 - 2 * wr * wi / 4);
        }
        int maxit = (int)(60 + 45 * log(zoom + 1));
        if (maxit > 5000) maxit = 5000;
        double logmax = log(maxit);

        #pragma omp parallel for schedule(dynamic, 1)
        for (int y = 0; y < H; y++) {
            char *p = rows + rowcap * y;
            int last = -1;
            double ci = CY + (y - H / 2.0) * scale * aspect;
            for (int x = 0; x < W; x++) {
                double cr = CX + (x - W / 2.0) * scale;
                /* supersample SS x SS points per character and average them;
                   inside points count as the densest character */
                double tsum = 0, nusum = 0;
                int nin = 0;
                for (int sy = 0; sy < SS; sy++)
                    for (int sx = 0; sx < SS; sx++) {
                        double nu;
                        if (escapes(cr + (sx + 0.5 - SS / 2.0) / SS * scale,
                                    ci + (sy + 0.5 - SS / 2.0) / SS * scale * aspect,
                                    maxit, &nu)) {
                            tsum += ramp_pos(nu, logmax);
                            nusum += nu;
                        } else {
                            tsum += 1;
                            nin++;
                        }
                    }
                char c;
                int col;
                if (nin * 4 >= SS * SS * 3) { /* mostly inside */
                    c = '@';
                    col = -2;
                    nubuf[y * W + x] = -1;
                } else {
                    double t = tsum / (SS * SS);
                    double nu = nusum / (SS * SS - nin);
                    tbuf[y * W + x] = (float)t;
                    c = ramp[(int)(t * (n - 1))];
                    col = (int)(nu * 2 + frame) % 64;
                    nubuf[y * W + x] = (float)nu;
                }
                if (COLOR && col != last) {
                    if (col == -2) { memcpy(p, inside_col, sizeof inside_col - 1); p += sizeof inside_col - 1; }
                    else { memcpy(p, pal[col], pal_len[col]); p += pal_len[col]; }
                    last = col;
                }
                *p++ = c;
            }
            *p++ = '\n';
            rowlen[y] = (int)(p - (rows + rowcap * y));
        }

        char *o = out;
        memcpy(o, "\033[H", 3); o += 3;
        for (int y = 0; y < H; y++) { memcpy(o, rows + rowcap * y, rowlen[y]); o += rowlen[y]; }
        char status[512], cyc[48] = "";
        if (steps) {
            if (rewinding) snprintf(cyc, sizeof cyc, " | rewinding");
            else snprintf(cyc, sizeof cyc, " | step %d/%d", step, steps);
        }
        snprintf(status, sizeof status,
                 " %s | zoom %.3gx | iter %d | %.0f fps%s%s%s | c = %.15f %+.15fi | q to quit",
                 label, zoom, maxit, fps, cyc, autopilot ? " | autopilot" : "",
                 zspeed < 0 ? " (backing out)" : zspeed < 0.99 ? " (steering)" : "",
                 ftype == F_JULIA ? jr : CX, ftype == F_JULIA ? ji : CY);
        o += sprintf(o, "\033[0m\033[7m%-*.*s\033[0m", W, W, status);
        fwrite(out, 1, o - out, stdout);
        fflush(stdout);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double ns = (t1.tv_sec - t0.tv_sec) * 1e9 + (t1.tv_nsec - t0.tv_nsec);
        if (ns < frame_ns) {
            struct timespec d = {0, (long)(frame_ns - ns)};
            nanosleep(&d, NULL);
            ns = frame_ns;
        }
        fps = 1e9 / ns;
        dt = ns / 1e9;
        if (dt > 0.2) dt = 0.2;

        /* Build the escape-time quantile table for the next frame */
        {
            long total = (long)W * H, stride = total / 8192 + 1;
            int ns = 0;
            for (long i = 0; i < total && ns < 8192; i += stride)
                if (nubuf[i] >= 0) samp[ns++] = nubuf[i];
            have_q = ns >= 32;
            if (have_q) {
                qsort(samp, ns, sizeof(float), cmpf);
                /* blend into the previous table so the character mapping
                   drifts smoothly instead of flickering frame to frame */
                double a = qready ? 1 - exp(-dt * 3) : 1;
                for (int k = 0; k <= NQ; k++)
                    qtab[k] += (samp[(long)k * (ns - 1) / NQ] - qtab[k]) * a;
                qready = 1;
            }
        }

        /* Autopilot: the interesting structure (spirals, mini-brots) lives on
           the edge of the set, so score each block of the screen by how
           close it is to ~20% black, plus how near to the set its outside
           points are, weighted towards the centre; steer to the best one.
           The zoom slows, stops or reverses when the screen gets too black. */
        zspeed = 1;
        if (autopilot) {
            const int BW = W >= 32 ? 16 : W / 2, BH = H >= 16 ? 8 : (H > 1 ? H / 2 : 1);
            int nbx = W / BW, nby = H / BH;
            double best = -1, bpx = W / 2.0, bpy = H / 2.0, cur = 0;
            double tx = (TR - CX) / scale + W / 2.0;
            double ty = (TI - CY) / (scale * aspect) + H / 2.0;
            for (int by = 0; by < nby; by++)
                for (int bx = 0; bx < nbx; bx++) {
                    double tsum = 0;
                    int ins = 0;
                    for (int y = by * BH; y < (by + 1) * BH; y++)
                        for (int x = bx * BW; x < (bx + 1) * BW; x++) {
                            if (nubuf[y * W + x] < 0) ins++;
                            else tsum += tbuf[y * W + x];
                        }
                    int cells = BW * BH;
                    double fin = (double)ins / cells;
                    double tmean = ins < cells ? tsum / (cells - ins) : 0;
                    /* sqrt(f)*(1-f)^2 peaks at f = 0.2; scaled to 1 there */
                    double edge = 3.5 * sqrt(fin) * (1 - fin) * (1 - fin);
                    double px = (bx + 0.5) * BW, py = (by + 0.5) * BH;
                    double dx = (px - W / 2.0) / (W / 2.0), dy = (py - H / 2.0) / (H / 2.0);
                    double sc = (edge + 0.3 * tmean) / (1 + 4 * (dx * dx + dy * dy));
                    if (sc > best) { best = sc; bpx = px; bpy = py; }
                    if (tx >= bx * BW && tx < (bx + 1) * BW && ty >= by * BH && ty < (by + 1) * BH)
                        cur = sc;
                }
            /* how black is the block-sized area we are zooming into? */
            long cin = 0, ccount = 0;
            for (int y = (H - BH) / 2; y < (H + BH) / 2; y++)
                for (int x = (W - BW) / 2; x < (W + BW) / 2; x++, ccount++)
                    cin += nubuf[y * W + x] < 0;
            black = ccount ? (double)cin / ccount : 0;
            /* hysteresis: only switch target when clearly better; a black
               centre makes switching easy */
            if (best > 2 * cur * (1 - black)) {
                TR = CX + (bpx - W / 2.0) * scale;
                TI = CY + (bpy - H / 2.0) * scale * aspect;
            }
            double k = 1 - exp(-dt * 1.2); /* smooth pan, ~0.8 s time constant */
            CX += (TR - CX) * k;
            CY += (TI - CY) * k;
            /* a smooth edge is ~50% black: fine. Mostly black = heading inside */
            if (black > 0.6) zspeed = black < 0.85 ? (0.85 - black) / 0.25 : -(black - 0.85) / 0.15;
        }
        /* ease the zoom speed: gentle start, and no sudden autopilot
           brakes, so the motion stays fluid */
        elapsed += dt;
        double ease = elapsed < 2 ? elapsed / 2 * (elapsed / 2) * (3 - elapsed) : 1;
        zsmooth += (zspeed * ease - zsmooth) * (1 - exp(-dt * 2));
        if (steps) {
            /* fixed zoom per step, so the rewind retraces the zoom-in;
               velocity is smoothed, giving soft turns at both ends */
            double vt;
            if (!rewinding) {
                vt = dz * zspeed;
                if (++step >= steps || (L >= Lmax && ftype != F_JULIA)) rewinding = 1;
            }
            if (rewinding) {
                double r = L / (60 * dz); /* slow down near the start */
                vt = -3 * dz * (r > 1 ? 1 : r < 0.05 ? 0.05 : r);
            }
            v += (vt - v) * 0.08;
            L += v;
            if (L > Lmax) L = Lmax;
            if (rewinding && L <= 0) break; /* back at the start: next pattern */
            if (L < 0) L = 0;
            scale = scale0 * pow(10, -L);
        } else {
            scale *= pow(10, -zsmooth * dt / zsecs);
            if (scale < scale0 * pow(10, -Lmax)) scale = scale0 * pow(10, -Lmax);
        }
        if (scale > scale0) scale = scale0;
        zoom = scale0 / scale;
        frame++;
    }
    return 0;
}
EOF
)

BIN=""
CC=$(command -v cc || command -v gcc || command -v clang)
if [[ -n $CC ]]; then
  CACHE="${XDG_CACHE_HOME:-$HOME/.cache}/mandelbrot"
  BIN="$CACHE/render-$(printf '%s' "$C_SRC" | md5sum | cut -c1-12)"
  if [[ ! -x $BIN ]]; then
    mkdir -p "$CACHE"
    printf '%s\n' "$C_SRC" > "$CACHE/render.c"
    { "$CC" -O3 -march=native -fopenmp -o "$BIN.tmp" "$CACHE/render.c" -lm 2>/dev/null ||
      "$CC" -O3 -o "$BIN.tmp" "$CACHE/render.c" -lm 2>/dev/null; } &&
      mv "$BIN.tmp" "$BIN" || BIN=""
  fi
fi

# A terminal's font can't be changed from inside it, so re-launch this script
# in a new fullscreen foot window with a small font.
if (( NEWWIN )) && [[ -z $MANDEL_CHILD && -n $WAYLAND_DISPLAY ]] && command -v foot >/dev/null; then
  MANDEL_CHILD=1 exec foot --log-level=error --fullscreen --app-id=mandelbrot \
    -o "font=JetBrainsMono Nerd Font:size=$FONTSIZE" "$(realpath "$0")" "${ARGS[@]}"
fi

# Zoom targets are Misiurewicz points (or their equivalent for the other
# fractals): exact points on the edge of the set where it is self-similar,
# so a straight zoom never ends up in the black inside and the patterns
# keep repeating forever. Format: fractal-type re im name
FNAMES=("" "Mandelbrot" "Multibrot z^3" "Multibrot z^4" "Burning Ship" "Tricorn" "Julia")
TARGETS=(
  "1 -0.77568376800905381 0.13646736829469011 Seahorse spiral"
  "1 0.29447552362463314 -0.01631644173420511 Elephant Valley"
  "1 -0.77468081302823577 0.13741680715479157 Double spiral"
  "1 -0.10109636384562215 0.95628651080914151 Quad spiral arms"
  "1 -1.54368901269207637 0 Antenna feather"
  "1 -0.10110375797210418 0.95629400795235409 Dendrite branch"
  "2 -0.54128766541271189 -0.64765741967952917 Triple spirals"
  "2 0.64565820918483896 0.29849474993225139 Three-way branch"
  "2 -0.32815751084797035 -0.73978523157516551 Trident"
  "3 -1.13729702786970988 0.27390253255809427 Four-fold spiral"
  "3 0.50511762029564145 1.01621834049617599 Cross branch"
  "3 -0.96261449128744214 0.22891784700162232 Quad bulb edge"
  "4 -1.57585866995559498 -0.03774438901365196 The armada"
  "4 -1.63178452920640549 -0.06282823000779364 Mast and rigging"
  "4 -1.86075577326614083 -0.01164947335319629 Antenna ships"
  "5 0.68637645986783169 0.87034206189472807 Broken shoreline"
  "5 -1.09692656551687873 0.15924841985780847 Wavy spiral"
  "5 0.53857134382105287 0.65482522817726196 Tricorn filaments"
  "6 0 0 Morphing"
)
# the targets belonging to the chosen fractal (FRACTAL 0 = all of them)
build_active() {
  ACTIVE=()
  local i f
  for i in "${!TARGETS[@]}"; do
    read -r f _ <<<"${TARGETS[$i]}"
    (( FRACTAL == 0 || f == FRACTAL )) && ACTIVE+=("$i")
  done
}
set_target() {
  local idx=${ACTIVE[$(( (TARGET - 1) % ${#ACTIVE[@]} ))]}
  read -r FTYPE CX CY TNAME <<<"${TARGETS[$idx]}"
  LABEL="${FNAMES[FTYPE]}: $TNAME"
}
build_active
if (( ${#ACTIVE[@]} == 0 )); then echo "unknown fractal -x $FRACTAL (use 0-6)" >&2; exit 1; fi
set_target

AWK=$(command -v mawk || command -v gawk || command -v awk)

cleanup() {
  kill "$RENDER_PID" 2>/dev/null
  printf '\033[0m\033[?25h\033[?1049l'
  stty "$OLD_STTY" 2>/dev/null
  exit 0
}
trap cleanup INT TERM EXIT
OLD_STTY=$(stty -g 2>/dev/null)
stty -echo -icanon min 0 time 0 2>/dev/null

# Alternate screen, hide cursor, clear
printf '\033[?1049h\033[?25l\033[2J'

render() {
  local rows cols
  read -r rows cols < <(stty size </dev/tty 2>/dev/null)
  [[ $cols =~ ^[0-9]+$ && $rows =~ ^[0-9]+$ ]] || { cols=$(tput cols); rows=$(tput lines); }
  [[ -n $BIN ]] && exec "$BIN" "$cols" "$rows" "$CX" "$CY" "$COLOR" "$ZSECS" "$AUTOPILOT" "$SS" "$RAMP" "$STEPS" "$FTYPE" "$LABEL"
  exec "$AWK" -v RAMP="$RAMP" -v ZSECS="$ZSECS" -v W="$cols" -v H="$rows" -v CX="$CX" -v CY="$CY" -v COLOR="$COLOR" '
  BEGIN {
    # Character ramp ordered from sparse (far from set) to dense (at the edge)
    ramp = " .\x27`^\",:;l!i><~+_-?][}{1)(|\\/tfjrxnuvcz0mwqpdbkhao*#&8%@$"
    if (RAMP != "") ramp = RAMP
    n = length(ramp)
    for (i = 1; i <= n; i++) ch[i-1] = substr(ramp, i, 1)
    inside = "@"
    H--                       # keep last line for the status bar
    aspect = 2.1              # terminal cells are ~2x taller than wide
    ln2 = log(2)

    # Precompute a smooth 256-colour palette (cube colours)
    for (i = 0; i < 64; i++) {
      t = i / 64
      r = int(5 * (0.5 + 0.5 * sin(6.2832 * (t + 0.00))))
      g = int(5 * (0.5 + 0.5 * sin(6.2832 * (t + 0.33))))
      b = int(5 * (0.5 + 0.5 * sin(6.2832 * (t + 0.67))))
      pal[i] = sprintf("\033[38;5;%dm", 16 + 36*r + 6*g + b)
    }

    scale = 3.5 / W           # complex units per column
    zf = exp(-log(10) / (ZSECS * 20))   # ~20 fps assumed for the awk renderer
    zoom = 1; frame = 0
    while (scale > 2e-15) {   # stop at double precision limit
      maxit = int(60 + 45 * log(zoom + 1))
      if (maxit > 3000) maxit = 3000
      logmax = log(maxit)
      out = "\033[H"
      for (y = 0; y < H; y++) {
        ci = CY + (y - H / 2) * scale * aspect
        line = ""; last = -1
        for (x = 0; x < W; x++) {
          cr = CX + (x - W / 2) * scale
          # quick reject: main cardioid and period-2 bulb
          q = (cr - 0.25) * (cr - 0.25) + ci * ci
          if (q * (q + cr - 0.25) < 0.25 * ci * ci || (cr + 1) * (cr + 1) + ci * ci < 0.0625) {
            it = maxit
          } else {
            zr = 0; zi = 0; zr2 = 0; zi2 = 0; it = 0
            while (zr2 + zi2 < 256 && it < maxit) {
              zi = 2 * zr * zi + ci
              zr = zr2 - zi2 + cr
              zr2 = zr * zr; zi2 = zi * zi
              it++
            }
          }
          if (it >= maxit) {
            c = inside; col = -2
          } else {
            # smooth (fractional) iteration count
            nu = it + 1 - log(log(zr2 + zi2) / 2 / ln2) / ln2
            if (nu < 1) nu = 1
            t = log(nu) / logmax
            if (t > 1) t = 1
            c = ch[int(t * (n - 1))]
            col = int(nu * 2 + frame) % 64
          }
          if (COLOR) {
            if (col != last) {
              line = line (col == -2 ? "\033[38;5;232m" : pal[col])
              last = col
            }
          }
          line = line c
        }
        out = out line "\n"
      }
      status = sprintf(" zoom %.3gx | iter %d | c = %.15f %+.15fi | q to quit", zoom, maxit, CX, CY)
      printf "%s\033[0m\033[7m%-" W "." W "s\033[0m", out, status
      fflush()
      if (maxit < 400) system("sleep 0.03")   # cap speed at shallow zoom
      scale *= zf; zoom /= zf; frame++
    }
  }'
}

# Restart the render at the new size whenever the window is resized
on_resize() { RESIZED=1; kill "$RENDER_PID" 2>/dev/null; }
trap on_resize WINCH

sleep 0.3   # let a fresh fullscreen window settle to its final size
while true; do
  printf '\033[0m\033[2J'
  render &
  RENDER_PID=$!
  RESIZED=0
  while kill -0 "$RENDER_PID" 2>/dev/null; do
    if read -rsn1 -t 0.2 key && [[ $key == q ]]; then cleanup; fi
  done
  # zoom finished: move on to the next target for some variety
  (( RESIZED )) || { TARGET=$(( TARGET % ${#ACTIVE[@]} + 1 )); set_target; }
done

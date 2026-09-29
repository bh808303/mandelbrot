# mandelbrot.sh

A full-screen, animated ASCII fractal zoom for the terminal, in a single bash script.

It zooms into hand-picked spots of the Mandelbrot set (Seahorse spiral, Elephant Valley, …) and five other fractals, with colour, smooth shading and a status bar showing what's playing.

## Features

- **Six fractals:** Mandelbrot, Multibrot z³ and z⁴, Burning Ship, Tricorn, and a morphing Julia set. You can also play all of them one after another.
- **Curated zoom targets:** six for the Mandelbrot set and three for each of the others.
- **Cycle mode (default):** zooms in, rewinds back out, then moves on to the next target.
- **One-long-zoom mode:** zooms down to the limit of floating-point precision.
- **Autopilot:** wanders around, steering towards detail and away from the black interior.
- **Anti-aliasing:** 1×1 up to 4×4 samples per character.
- **Fast renderer:** a small C program is embedded in the script and compiled on first run. It uses all CPU cores through OpenMP. If there's no C compiler, the script falls back to a pure awk renderer.

## Requirements

- bash and a terminal that supports ANSI colours
- Optional: `cc`, `gcc` or `clang`, for the fast renderer
- Optional: [foot](https://codeberg.org/dnkl/foot) on Wayland, for the `-w` fullscreen small-font window

## Usage

```sh
git clone https://github.com/bh808303/mandelbrot.git
cd mandelbrot
./mandelbrot.sh
```

Quit with `q` or `Ctrl+C`.

```
./mandelbrot.sh [-x fractal] [-m] [-t target] [-z secs] [-a] [-c|-o] [-l steps] [-s n] [-r] [-w] [-f size]
```

| Option     | Description |
|------------|-------------|
| `-x N`     | Fractal: `1` Mandelbrot (default), `2` Multibrot z³, `3` Multibrot z⁴, `4` Burning Ship, `5` Tricorn, `6` Julia (morphing), `0` all of them in turn |
| `-m`       | Monochrome (no ANSI colours) |
| `-t N`     | Start at the fractal's Nth zoom target (default `1`) |
| `-z SECS`  | Zoom speed in seconds per 10× zoom (default `6`; higher is slower) |
| `-a`       | Autopilot: wander towards detail instead of zooming straight in |
| `-c`       | Cycle (default): zoom in for a fixed number of steps, rewind, then go to the next target |
| `-o`       | One long zoom per target, down to the precision limit |
| `-l STEPS` | Steps (frames) per zoom-in when cycling (default `1000`) |
| `-s N`     | Smoothing: N×N samples per character, `1`–`4` (default `3`) |
| `-r`       | Use the long 57-character ramp instead of the default ` .,:-~=+*o#%&8@` |
| `-w`       | Open a new fullscreen foot window with a small font |
| `-f SIZE`  | Font size for `-w` (default `5`) |

Mandelbrot zoom targets for `-t`: 1 Seahorse spiral, 2 Elephant Valley, 3 Double spiral, 4 Quad spiral arms, 5 Antenna feather, 6 Dendrite branch.

### Examples

```sh
./mandelbrot.sh -w                # high-res fullscreen window (foot)
./mandelbrot.sh -x 0 -a           # every fractal, on autopilot
./mandelbrot.sh -x 4 -o -z 10     # one slow, deep zoom into the Burning Ship
./mandelbrot.sh -s 1 -m           # fastest: no smoothing, no colour
```

## Notes

The compiled renderer is cached in `${XDG_CACHE_HOME:-~/.cache}/mandelbrot/`. It's rebuilt automatically when the embedded C source changes.

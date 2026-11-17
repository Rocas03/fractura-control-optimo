#!/bin/sh
# Regenera la ENTREGA FINAL: los json de final/resultados/, las figuras de final/figuras/, y
# comprueba que cada figura compile SOLA. Uso, desde cualquier carpeta: sh final/regenerar.sh
# `avance/` está entregado (26-oct) y no se toca: lo regenera avance/regenerar.sh, y que siga dando
# los mismos bytes es el control de que no rompimos lo entregado.
# Julia NO se busca en el PATH: el envoltorio julia/bin/jl fija 1.11.4 (el PATH trae 1.13, que al
# resolver el Manifest sale a la red y el sandbox la corta).
set -eu
cd "$(dirname "$0")/.."

for herramienta in latexmk xelatex pdfinfo; do
    command -v "$herramienta" >/dev/null 2>&1 || { echo "final/regenerar.sh: falta $herramienta" >&2; exit 2; }
done

echo "1/2 paquete Fractura: json y figuras de final/"
julia/bin/jl --project=julia julia/regenerar_final.jl

echo "2/2 cada figura, compilada sola"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/final-figuras.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
cat > "$TMP/marco.tex" <<'TEX'
\documentclass[10pt,a4paper]{article}
\usepackage{fontspec}
\usepackage[spanish,es-nodecimaldot,es-noshorthands]{babel}
\usepackage[a4paper,margin=1.6cm]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{pgfplots}\pgfplotsset{compat=1.18}\usepgfplotslibrary{groupplots}
\tracinglostchars=3
\pagestyle{empty}
\begin{document}
\noindent\input{FIGURA}
\end{document}
TEX

fallas=0
for f in final/figuras/*.tex; do
    nombre=$(basename "$f" .tex)
    cp "$f" "$TMP/$nombre.tex"
    sed "s/FIGURA/$nombre/" "$TMP/marco.tex" > "$TMP/marco-$nombre.tex"
    if ! (cd "$TMP" && latexmk -xelatex -g -interaction=nonstopmode -halt-on-error \
            "marco-$nombre.tex" >"$TMP/$nombre.salida" 2>&1); then
        tail -20 "$TMP/$nombre.salida" >&2
        echo "final/regenerar.sh: $nombre.tex NO compila" >&2
        fallas=$((fallas + 1))
        continue
    fi
    sin_glifo=$(grep -c 'Missing character' "$TMP/marco-$nombre.log" || true)
    cajas=$(grep -c -E 'Overfull|Underfull' "$TMP/marco-$nombre.log" || true)
    paginas=$(pdfinfo "$TMP/marco-$nombre.pdf" | awk '/^Pages:/ {print $2}')
    echo "  $nombre.tex: rc=0, $paginas pág., $sin_glifo caracteres sin glifo, $cajas cajas Overfull/Underfull"
    if [ "$sin_glifo" -ne 0 ] || [ "$cajas" -ne 0 ]; then
        fallas=$((fallas + 1))
    fi
done

if [ "$fallas" -ne 0 ]; then
    echo "final/regenerar.sh: $fallas figura(s) no cumplen (0 sin glifo, 0 cajas)" >&2
    exit 1
fi
echo "listo: final/resultados/ y final/figuras/"

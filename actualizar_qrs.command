#!/bin/bash
cd "$(dirname "$0")"
python3 generar_qr.py
echo ""
read -p "Presiona ENTER para cerrar esta ventana..."

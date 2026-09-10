import sys

EXP = [0] * 512
LOG = [0] * 256
x = 1
for i in range(255):
    EXP[i] = x
    EXP[i + 255] = x
    LOG[x] = i
    x = (x << 1) ^ (0x11d if (x & 0x80) else 0)

def gmult(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[(LOG[a] + LOG[b]) % 255]

def rs_generator_poly(nsym):
    g = [1]
    for i in range(nsym):
        ng = [0] * (len(g) + 1)
        root = EXP[i]
        for j, c in enumerate(g):
            ng[j] ^= gmult(c, root)
            ng[j + 1] ^= c
        g = ng
    return g

def rs_encode(msg, nsym):
    gen = rs_generator_poly(nsym)
    msg_out = list(msg) + [0] * nsym
    for i in range(len(msg)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(1, len(gen)):
                msg_out[i + j] ^= gmult(gen[len(gen) - 1 - j], coef)
    return msg_out[len(msg):]

TABLE_L = {
    1: (21, 26, [(19, 7, 1)], []),
    2: (25, 44, [(34, 10, 1)], [6, 18]),
    3: (29, 70, [(55, 15, 1)], [6, 22]),
    4: (33, 100, [(80, 20, 1)], [6, 26]),
    5: (37, 134, [(108, 26, 1)], [6, 30]),
    6: (41, 172, [(68, 18, 2)], [6, 34]),
}

TABLE_M = {
    1: (21, 26, [(16, 10, 1)], []),
    2: (25, 44, [(28, 16, 1)], [6, 18]),
    3: (29, 70, [(44, 26, 1)], [6, 22]),
    4: (33, 100, [(32, 18, 2)], [6, 26]),
    5: (37, 134, [(43, 24, 2)], [6, 30]),
    6: (41, 172, [(27, 16, 4)], [6, 34]),
}

def get_format_bits(ec_level_bits, mask):
    val = (ec_level_bits << 3) | mask
    rem = val << 10
    for i in range(14, 9, -1):
        if rem & (1 << i):
            rem ^= (0x537 << (i - 10))
    bits = (val << 10) | rem
    return bits ^ 0x5412

def create_qr(text, ec_level='M'):
    raw_data = text.encode('utf-8')
    data_len = len(raw_data)
    
    table = TABLE_M if ec_level == 'M' else TABLE_L
    ec_bits = 0b00 if ec_level == 'M' else 0b01
    
    version = None
    for v in range(1, 7):
        size, total_cw, blocks, aligns = table[v]
        total_data_cw = sum(data_cw * count for data_cw, ec_cw, count in blocks)
        if data_len + 3 <= total_data_cw:
            version = v
            break
            
    if version is None:
        raise ValueError("Data too long for versions 1-6")
        
    size, total_cw, blocks, aligns = table[version]
    total_data_cw = sum(data_cw * count for data_cw, ec_cw, count in blocks)
    
    # Encode bit stream
    # Mode: 0100 (Byte)
    bits = "0100"
    # Char count indicator (8 bits for v1-9)
    bits += f"{data_len:08b}"
    for b in raw_data:
        bits += f"{b:08b}"
    # Terminator
    bits += "0000"[:min(4, total_data_cw * 8 - len(bits))]
    # Pad to multiple of 8
    if len(bits) % 8 != 0:
        bits += "0" * (8 - (len(bits) % 8))
    # Pad bytes
    data_bytes = bytearray()
    for i in range(0, len(bits), 8):
        data_bytes.append(int(bits[i:i+8], 2))
        
    pad_toggle = False
    while len(data_bytes) < total_data_cw:
        data_bytes.append(0x11 if pad_toggle else 0xEC)
        pad_toggle = not pad_toggle
        
    # Block splitting & Reed-Solomon
    block_data = []
    block_ec = []
    offset = 0
    for data_cw, ec_cw, count in blocks:
        for _ in range(count):
            chunk = data_bytes[offset:offset+data_cw]
            offset += data_cw
            ec = rs_encode(chunk, ec_cw)
            block_data.append(chunk)
            block_ec.append(ec)
            
    # Interleave data
    interleaved = bytearray()
    max_dlen = max(len(b) for b in block_data)
    for i in range(max_dlen):
        for b in block_data:
            if i < len(b):
                interleaved.append(b[i])
    max_eclen = max(len(e) for e in block_ec)
    for i in range(max_eclen):
        for e in block_ec:
            if i < len(e):
                interleaved.append(e[i])
                
    final_bits = ""
    for b in interleaved:
        final_bits += f"{b:08b}"
        
    # Remainder bits (for v2-6)
    rem_bits_count = {1: 0, 2: 7, 3: 7, 4: 7, 5: 7, 6: 7}
    final_bits += "0" * rem_bits_count[version]
    
    # Matrix setup
    matrix = [[None] * size for _ in range(size)]
    reserved = [[False] * size for _ in range(size)]
    
    # 1. Finder patterns
    def set_finder(orow, ocol):
        for r in range(-1, 8):
            for c in range(-1, 8):
                nr, nc = orow + r, ocol + c
                if 0 <= nr < size and 0 <= nc < size:
                    reserved[nr][nc] = True
                    if 0 <= r <= 6 and 0 <= c <= 6:
                        if r in (0, 6) or c in (0, 6) or (2 <= r <= 4 and 2 <= c <= 4):
                            matrix[nr][nc] = 1
                        else:
                            matrix[nr][nc] = 0
                    else:
                        matrix[nr][nc] = 0
                        
    set_finder(0, 0)
    set_finder(0, size - 7)
    set_finder(size - 7, 0)
    
    # 2. Timing patterns
    for i in range(size):
        if matrix[6][i] is None:
            matrix[6][i] = 1 if i % 2 == 0 else 0
            reserved[6][i] = True
        if matrix[i][6] is None:
            matrix[i][6] = 1 if i % 2 == 0 else 0
            reserved[i][6] = True
            
    # 3. Alignment patterns
    if len(aligns) > 0:
        for ar in aligns:
            for ac in aligns:
                # Do not place on finders
                if (ar < 9 and ac < 9) or (ar < 9 and ac >= size - 9) or (ar >= size - 9 and ac < 9):
                    continue
                for r in range(-2, 3):
                    for c in range(-2, 3):
                        nr, nc = ar + r, ac + c
                        reserved[nr][nc] = True
                        if r in (-2, 2) or c in (-2, 2) or (r == 0 and c == 0):
                            matrix[nr][nc] = 1
                        else:
                            matrix[nr][nc] = 0
                            
    # 4. Dark module
    matrix[size - 8][8] = 1
    reserved[size - 8][8] = True
    
    # 5. Format info positions (reserve)
    for i in range(9):
        reserved[8][i] = True
        reserved[i][8] = True
    for i in range(8):
        reserved[8][size - 1 - i] = True
        reserved[size - 1 - i][8] = True
        
    # Mask functions
    mask_funcs = [
        lambda r, c: (r + c) % 2 == 0,
        lambda r, c: r % 2 == 0,
        lambda r, c: c % 3 == 0,
        lambda r, c: (r + c) % 3 == 0,
        lambda r, c: (r // 2 + c // 3) % 2 == 0,
        lambda r, c: ((r * c) % 2) + ((r * c) % 3) == 0,
        lambda r, c: (((r * c) % 2) + ((r * c) % 3)) % 2 == 0,
        lambda r, c: (((r + c) % 2) + ((r * c) % 3)) % 2 == 0,
    ]
    
    best_mask = 0
    best_score = float('inf')
    best_grid = None
    
    for mask_idx in range(8):
        grid = [row[:] for row in matrix]
        bit_idx = 0
        
        # Place data
        col = size - 1
        up = True
        while col > 0:
            if col == 6:
                col -= 1
            row_range = range(size - 1, -1, -1) if up else range(size)
            for r in row_range:
                for c in (col, col - 1):
                    if not reserved[r][c]:
                        bit_val = int(final_bits[bit_idx]) if bit_idx < len(final_bits) else 0
                        bit_idx += 1
                        if mask_funcs[mask_idx](r, c):
                            bit_val ^= 1
                        grid[r][c] = bit_val
            up = not up
            col -= 2
            
        # Place format bits
        fbits = get_format_bits(ec_bits, mask_idx)
        # Around top-left
        coords_tl = [
            (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5),
            (8, 7), (8, 8), (7, 8), (5, 8), (4, 8), (3, 8), (2, 8), (1, 8), (0, 8)
        ]
        for i, (r, c) in enumerate(coords_tl):
            grid[r][c] = (fbits >> (14 - i)) & 1
            
        # Top-right and bottom-left
        for i in range(7):
            grid[size - 1 - i][8] = (fbits >> (14 - i)) & 1
        for i in range(8):
            grid[8][size - 8 + i] = (fbits >> (7 - i)) & 1
            
        # Evaluate penalty
        score = 0
        # N1: 5 in a row
        for r in range(size):
            run_color, run_len = None, 0
            for c in range(size):
                v = grid[r][c]
                if v == run_color:
                    run_len += 1
                else:
                    if run_len >= 5:
                        score += 3 + (run_len - 5)
                    run_color, run_len = v, 1
            if run_len >= 5:
                score += 3 + (run_len - 5)
                
        for c in range(size):
            run_color, run_len = None, 0
            for r in range(size):
                v = grid[r][c]
                if v == run_color:
                    run_len += 1
                else:
                    if run_len >= 5:
                        score += 3 + (run_len - 5)
                    run_color, run_len = v, 1
            if run_len >= 5:
                score += 3 + (run_len - 5)
                
        # N2: 2x2 blocks
        for r in range(size - 1):
            for c in range(size - 1):
                if grid[r][c] == grid[r+1][c] == grid[r][c+1] == grid[r+1][c+1]:
                    score += 3
                    
        # N4: Dark ratio
        dark_count = sum(sum(row) for row in grid)
        pct = (dark_count * 100) // (size * size)
        score += abs(pct - 50) // 5 * 10
        
        if score < best_score:
            best_score = score
            best_mask = mask_idx
            best_grid = grid
            
    return best_grid

def qr_to_svg(grid, border=4, scale=10):
    size = len(grid)
    total_size = (size + border * 2) * scale
    
    rects = []
    for r in range(size):
        for c in range(size):
            if grid[r][c] == 1:
                x = (c + border) * scale
                y = (r + border) * scale
                rects.append(f'<rect x="{x}" y="{y}" width="{scale}" height="{scale}" fill="#1e293b"/>')
                
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_size} {total_size}" width="{total_size}" height="{total_size}">
  <rect width="{total_size}" height="{total_size}" fill="#ffffff" rx="16"/>
  {''.join(rects)}
</svg>'''
    return svg

def grid_to_png(grid, border=4, scale=14, fg=(15, 23, 42), bg=(255, 255, 255)):
    import zlib
    import struct
    size = len(grid)
    width = (size + border * 2) * scale
    height = width
    raw_data = bytearray()
    for y in range(height):
        raw_data.append(0)
        grid_y = y // scale - border
        for x in range(width):
            grid_x = x // scale - border
            if 0 <= grid_y < size and 0 <= grid_x < size and grid[grid_y][grid_x] == 1:
                raw_data.extend(fg)
            else:
                raw_data.extend(bg)
    png = bytearray(b'\x89PNG\r\n\x1a\n')
    def write_chunk(ctype, data):
        png.extend(struct.pack('>I', len(data)))
        png.extend(ctype)
        png.extend(data)
        crc = zlib.crc32(ctype + data) & 0xffffffff
        png.extend(struct.pack('>I', crc))
    write_chunk(b'IHDR', struct.pack('>IIBBBBB', width, height, 8, 2, 0, 0, 0))
    write_chunk(b'IDAT', zlib.compress(raw_data, 6))
    write_chunk(b'IEND', b'')
    return bytes(png)

# ============================================================
# CONFIGURACIÓN: URL BASE DE TU PROYECTO PUBLICADO
# Puedes cambiar este enlace aquí mismo, o escribirlo cuando
# ejecutes este programa.
# ============================================================
URL_BASE_POR_DEFECTO = "https://ignis-registros.vercel.app"

if __name__ == '__main__':
    import os
    
    # Si se pasó la URL por la terminal (ej: python3 generar_qr.py "https://...")
    if len(sys.argv) > 1:
        base_url = sys.argv[1].strip()
    else:
        # Modo interactivo amigable
        if sys.stdin.isatty():
            print("\n" + "="*60)
            print("        ACTUALIZADOR DE CÓDIGOS QR DEL EQUIPO")
            print("="*60)
            print(f"URL actual predeterminada: {URL_BASE_POR_DEFECTO}")
            try:
                entrada = input("Enlace web: ").strip()
                base_url = entrada if entrada else URL_BASE_POR_DEFECTO
            except (EOFError, KeyboardInterrupt):
                base_url = URL_BASE_POR_DEFECTO
        else:
            base_url = URL_BASE_POR_DEFECTO

    base_url = base_url.rstrip("/")
    os.makedirs("qrs", exist_ok=True)
    
    integrantes = ["denisse", "mildred", "axel", "jacobo", "diego", "ximena", "dereck", "alexa", "jorge", "santiago"]
    print(f"\nGenerando códigos QR (SVG y PNG) con la ruta: {base_url}\n")
    
    for nombre in integrantes:
        url = f"{base_url}/{nombre}.html"
        grid = create_qr(url, 'M')
        
        # Guardar SVG
        svg = qr_to_svg(grid, border=4, scale=12)
        svg_filepath = os.path.join("qrs", f"qr-{nombre}.svg")
        with open(svg_filepath, "w", encoding="utf-8") as f:
            f.write(svg)
            
        # Guardar PNG (alta resolución, ideal para celulares e impresión)
        png_bytes = grid_to_png(grid, border=4, scale=14)
        png_filepath = os.path.join("qrs", f"qr-{nombre}.png")
        with open(png_filepath, "wb") as f:
            f.write(png_bytes)
            
        print(f"  [OK] SVG y PNG: {nombre} -> {url}")
        
    print("\n¡Listo! Todos los códigos QR en SVG y PNG fueron actualizados exitosamente en la carpeta qrs/.")
    print("Puedes abrir el archivo index.html para verlos o descargarlos.\n")



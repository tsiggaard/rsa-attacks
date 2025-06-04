import rsa
import random
import json

s_prev = 0
oracle_counter = 0
e = 65537

def floordiv(a, b):
  return a // b

def ceildiv(a, b):
  return -(-a // b)

def generate_key_pair(k):
  rsa_mod_size = k*8
  print("Working on RSA(", rsa_mod_size, ")")
  return rsa.newkeys(rsa_mod_size)

def generate_padded_message(m, k):
  # 00 02 padding 00 message
  # padding must be at least 8 bytes
  if k < 12:
    raise Exception("k must be at least 12 bytes")
  if m.bit_length() > (k - 11) * 8:
    raise Exception("Message too long or k too short")
  m_length = ceildiv(m.bit_length(), 8)
  padding_length = k - 3 - m_length
  padding = bytes([random.randint(1, 255) for _ in range(padding_length)])
  return int.from_bytes(b'\x00\x02' + padding + b'\x00' + m.to_bytes(m_length, byteorder='big'), byteorder='big')

def get_message_from_padded(m, k):
  m_b = m.to_bytes(k, byteorder='big')
  try:
    padding_end = m_b.index(0x00, 2)
    return str(m_b[padding_end + 1:], 'utf-8')
  except ValueError:
    return ""

def padding_check(m, k):
  global oracle_counter
  oracle_counter += 1
  m_b = m.to_bytes(k, byteorder='big')
  if (m_b[0] != 0x00) or (m_b[1] != 0x02):
    return False
  if k > 11:
    if any(x == 0x00 for x in m_b[2:10]):
      return False
  try:
    padding_end = m_b.index(0x00, 2)
  except ValueError:
    return False
  return padding_end < k - 1

def oracle(c, d, n, k):
  if d == 0:
    return padding_check(c, k)
  m = pow(c, d, n)
  return padding_check(m, k)

def union_range(M):
  b = []
  for begin,end in sorted(M):
    if b and b[-1][1] >= begin - 1:
      b[-1][1] = max(b[-1][1], end)
    else:
      b.append([begin, end])
  return b

def step_1(m, n, d, k, enc=False, debug=False):
  global s_prev
  if debug:
    print(" STEP 1 ")
  for s in range(1, n):
    se = pow(s, e, n)
    m_i = (m * se) % n if enc else (m * s) % n
    if oracle(m_i, d, n, k):
      if debug:
        print("s:", format(s, ","))
        print("Calls to oracle:", format(oracle_counter, ","))
        print("Padding conforming m:" + hex(m_i))
      s_prev = s
      return m_i, s

def step_2_a_b(m, n, d, k, n_3B, enc=False, debug=False):
  global s_prev
  if debug:
    print(" STEP 2 A B ")
  for s_i in range(max(s_prev + 1, n_3B), n):
    se = pow(s_i, e, n)
    m_i = (m * se) % n if enc else (m * s_i) % n
    if oracle(m_i, d, n, k):
      if debug:
        print("s_i:", format(s_i, ","))
        print("Calls to oracle:", format(oracle_counter, ","))
        print("Padding conforming m:" + hex(m_i))
      s_prev = s_i
      return

def step_2_c(m, n, d, k, B, M_prev, enc=False, debug=False):
  global s_prev
  if debug:
    print(" STEP 2 C ")
  found = False
  r = ceildiv(2 * (M_prev[0][1] * s_prev - B), n)
  while not found:
    m_min_r = ceildiv(2*B + r*n, M_prev[0][1])
    m_max_r = floordiv(3*B + r*n, M_prev[0][0]) + 1
    # if debug:
    #   print("Range:", m_min_r, "-", m_max_r)
    for s in range(m_min_r, m_max_r):
      se = pow(s, e, n)
      m_i = (m * se) % n if enc else (m * s) % n
      if oracle(m_i, d, n, k):
        if debug:
          print("r:", r)
          print("s_i:", format(s, ","))
          print("Padding conforming m:" + hex(m_i))
          print("M_prev:", M_prev)
        s_prev = s
        found = True
        break
    r += 1

def step_3(n, B, M_prev, debug=False):
  global s_prev
  if debug:
    print(" STEP 3 ")
  M_current = []
  for a, b in M_prev:
    r_min = ceildiv((a * s_prev - 3*B + 1), n)
    r_max = floordiv((b * s_prev - 2*B), n) + 1
    for r in range(r_min, r_max):
      m_min_r = ceildiv((2*B + r*n), s_prev)
      m_max_r = floordiv((3*B - 1 + r*n), s_prev)
      m_min = max(a, m_min_r)
      m_max = min(b, m_max_r)
      if m_min <= m_max:
        M_current.append([m_min, m_max])
  result = union_range(M_current)
  return result

def step_4(a, s0, n, debug=False):
  if debug:
    print(" STEP 4 ")
  s0_1 = pow(s0, -1, n)
  m = (a * s0_1) % n
  return m

def write_intervals_to_json(iteration, intervals, filename="intervals.json"):
    """Write the current iteration and its intervals to a JSON file."""
    data_entry = {iteration: [{"a": a, "b": b} for a, b in intervals]}
    try:
        with open(filename, "r") as f:
            data = json.load(f) 
    except (FileNotFoundError, json.JSONDecodeError):
        data = [] 

    data.append(data_entry)
    
    with open(filename, "w") as f:
      json.dump(data, f, indent=2)

def loop(m, k, n, d, enc=False, debug=False):
  global s_prev
  global oracle_counter
  step_2_a_b_counter = 0
  i = 1
  s_prev = 0
  oracle_counter = 0
  B = 2**(8*(k-2))
  n_3B = ceildiv(n, 3*B)
  M = [[2*B, 3*B-1]]
  m, s0 = step_1(m, n, d, k, enc, debug)
  if debug:
    print("equal?", M[0][0] != M[0][1])
    write_intervals_to_json(i, M)
  while not (len(M) == 1 and (M[0][0] == M[0][1])):
    if debug:
      print("LOOP, i:", i)
    if len(M) != 1 or i == 1:
      step_2_a_b(m, n, d, k, n_3B, enc, debug)
      step_2_a_b_counter += 1
    else:
      step_2_c(m, n, d, k, B, M, enc, debug)
    M = step_3(n, B, M, debug)
    i += 1
    if debug:
      write_intervals_to_json(i, M)
  return step_4(M[0][0], s0, n, debug), step_2_a_b_counter

def bleichenbacher(m, k, n=0, d=0, enc=False, debug=False):
  pubKey, privKey = generate_key_pair(k)
  n = pubKey.n if n == 0 else n
  if enc and d == 0:
    d = privKey.d
  elif not enc and d != 0:
    d = 0
  if debug:
    print("n:", n)
    print("d:", d)
  c = pow(m, e, n) if enc else m 
  res, step_2_a_b_counter = loop(c, k, n, d, enc, debug)
  if debug:
    print("s:", s_prev)
    print("Result: " + hex(res))
    print("Is start and end message the same?", m == res)
    print("Calls to oracle:", format(oracle_counter, ","))
    print("Calculated message:", get_message_from_padded(res, k))
  print("Finished with RSA(", k * 8, ")")
  return (res, oracle_counter, step_2_a_b_counter)

def example():
  print("\nStart example")
  m = "Johan Enevoldsen, Emil Fra Lønneberg og Thomas Siggaard"
  k=128
  m_padded = generate_padded_message(int.from_bytes(m.encode(), byteorder='big'), k)
  res = bleichenbacher(m_padded, k, debug=True)

example()
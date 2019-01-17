let doc = """
dump_mmmer

Usage:
  dump_mmmer [options] <fasta_file_name>

Options:
  -h --help                       Show this screen
  -w --windowsize <windowsize>    Window size [default: 64]
"""

import streams
import strfmt
import tables
import strutils
import sequtils
import docopt
import parseutils

let args = docopt(doc, version = "dump minimizer")

type
  pos_kmer = tuple[pos:uint32, kmer:uint32]

var 
  fn = $args["<fasta_file_name>"]
  ws = parseInt($args["--windowsize"]).uint32

echo fn

var
  fs = newFileStream(fn, fmRead)
  line = ""

# not used
#[
  base_to_code = {'A':0.uint32, 'C':1.uint32,
                  'G':2.uint32, 'T':3.uint32,
                  'a':0.uint32, 'c':1.uint32,
                  'g':2.uint32, 't':3.uint32,
                  'N':0.uint32}.toTable
]#

  code_to_base = {0.uint32:'A', 1.uint32:'C',
                  2.uint32:'G', 3.uint32:'T'}.toTable


let xor_key = 0x7ed55d16.uint32
let ksize = 16.uint32

var mask = 0xFFFFFFFF.uint32 shr (32.uint32 - ksize * 2)


proc rc_DNA_seq(dna_seq:var string) : void {.inline.} =
  var rc_map = {'A':'T', 'C':'G',
                'G':'C', 'T':'A',
                'a':'t', 'c':'g',
                'g':'c', 't':'a',
                'N':'N'}.toTable
  for i in 0..toInt(dna_seq.len/2-1):
    swap dna_seq[i], dna_seq[^(i+1)]
  for i in 0..<dna_seq.len:
    dna_seq[i] = rc_map[dna_seq[i]]


proc decode_hash(hashcode:uint32) : string  {.inline.} =
  var hc: uint32
  hc = hashcode xor xor_key
  hc = hc and mask
  var t: uint32
  var s = newSeq[char](0)
  for i in 0..<ksize:
    s.add(code_to_base[hc and 0x3])
    t = hc and 0x3
    #echo hashcode, " ", hashcode xor xor_key, " ",hc, " ",t
    hc = hc shr 2
  for i in 0..toInt(s.len/2-1):
    swap s[i], s[^(i+1)]
  return s.mapIt(string, $it).join


proc find_minimizers(dna_seq:string): void =
  var 
    c: char
    p_mer: uint32
    h_mer: uint32
    w_start: uint32
    w_end: uint32
    mmer_seq: seq[uint32]
    pos: uint32
    pos2: uint32
    pm: pos_kmer 
    c_minimizer : pos_kmer

  mmer_seq= newSeq[uint32](dna_seq.len)
  p_mer = 0x00000000.uint32
  pos = 0.uint32

  for c in dna_seq:
    p_mer = p_mer shl 2 
    case c
    of 'C', 'c':
      p_mer = p_mer or 1.uint32
    of 'G', 'g':
      p_mer = p_mer or 2.uint32
    of 'T', 't':
      p_mer = p_mer or 3.uint32
    else:
      p_mer = p_mer or 0.uint32

    p_mer = p_mer and mask
    h_mer = p_mer xor xor_key
    if pos.int - ksize.int >= 0:
      mmer_seq[pos.int - ksize.int] = h_mer
    inc(pos)

  c_minimizer.pos = 0.uint32
  c_minimizer.kmer = 0xFFFFFFFF.uint32 
  w_start = 0
  w_end = ws
  for pos in w_start ..< w_end:
    if mmer_seq[pos.int] < c_minimizer.kmer:
      c_minimizer.pos = pos
      c_minimizer.kmer = mmer_seq[pos.int]
  echo "0 ", c_minimizer.pos,  " ", c_minimizer.kmer, " ", decode_hash(c_minimizer.kmer)

  for pos in ws ..< dna_seq.len.uint32 - ksize:
    # echo "X ", pos.int, " ", mmer_seq[pos.int], " ", c_minimizer.kmer
    if mmer_seq[pos.int] < c_minimizer.kmer:
      c_minimizer.pos = pos
      c_minimizer.kmer = mmer_seq[pos.int]
      echo "1 ",c_minimizer.pos,  " ", c_minimizer.kmer, " ", decode_hash(c_minimizer.kmer)
      continue

    if pos.int - c_minimizer.pos.int >= ws.int:
      w_start = c_minimizer.pos + 1
      w_end = w_start + ws
      c_minimizer.kmer = 0xFFFFFFFF.uint32 
      for pos2 in w_start ..< w_end:
        if mmer_seq[pos2.int] < c_minimizer.kmer:
          c_minimizer.pos = pos2
          c_minimizer.kmer = mmer_seq[pos2.int]
      echo  "2 ", c_minimizer.pos,  " ", c_minimizer.kmer, " ", decode_hash(c_minimizer.kmer)


var 
  dna_seq: string
  seq_name: string


if not isNil(fs):
  while fs.readLine(line):
    if line[0] == '>':
      if not isNil(seq_name):
        if dna_seq.len < ws.int:
          seq_name = line.strip
          dna_seq = ""
          continue
        echo seq_name, "|", "n"
        find_minimizers(dna_seq)
        rc_DNA_seq(dna_seq)
        echo seq_name, "|", "c"
        find_minimizers(dna_seq)
      seq_name = line.strip
      dna_seq = ""
      continue

    if line[0] != '>':
      dna_seq.add(line.strip)

  echo seq_name, "|", "n"
  find_minimizers(dna_seq)
  rc_DNA_seq(dna_seq)
  echo seq_name, "|", "r"
  find_minimizers(dna_seq)
  fs.close()

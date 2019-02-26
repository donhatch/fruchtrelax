#!/usr/bin/python3

import math
import re
import sys
import time

def vxs(v,s): return [x*s for x in v]  # vector times scalar
def vpv(v,w): return [x+y for x,y in zip(v,w)]  # vector plus vector
def vmv(v,w): return [x-y for x,y in zip(v,w)]  # vector minus vector
def lerp(v0,v1,t): return [x0*(1-t)+x1*t for x0,x1 in zip(v0,v1)]  # linear interpolation
def dot(v,w): return sum(x*y for x,y in zip(v,w))  # vector dot product
def length2(v): return dot(v,v)  # vector length (norm) squared
def length(v): return math.sqrt(length2(v))  # vector length (norm)
def dist(v,w): return length(vmv(w,v))

def relaxOnAnySizeSphere(verts,
                         edgeLengthConstraints,
                         nIters):
  verboseLevel = 1  # set to 2 to show one '.' per iteration
  if verboseLevel >= 1: sys.stderr.write("    in relaxOnAnySizeSphere\n")
  if verboseLevel >= 1: sys.stderr.write("      initial verts = %r\n"%(verts,))
  if verboseLevel >= 1: sys.stderr.write("      edgeLengthConstraints = %r\n"%(edgeLengthConstraints,))
  if verboseLevel >= 1: sys.stderr.write("      nIters = %r\n"%(nIters,))
  if verboseLevel >= 1: sys.stderr.write("      iterating... ")
  if verboseLevel >= 1: sys.stderr.flush()

  t0 = time.time()
  for iIter in range(nIters):
    if verboseLevel >= 2: sys.stderr.write(".")
    if verboseLevel >= 2: sys.stderr.flush()
    # Perturb verts to try to approach edge length constraints.
    # If a vertex is being pushed/pulled by multiple edges,
    # move it to the average of where they want it to go.
    newVerts = [(0.,0.,0.) for vert in verts]
    nContributions = [0] * len(newVerts)
    for v0,v1,relation,targetEdgeLength in edgeLengthConstraints:
      currentEdgeLength = dist(verts[v0],verts[v1])
      if ((relation=="<=" and currentEdgeLength<targetEdgeLength) or
          (relation==">=" and currentEdgeLength>targetEdgeLength)):
        # this spring already satisfies the constraint and it still will if perturbed a little;
        # so don't contribute anything.
        continue
      # this spring wants to change (or keep) its length from currentEdgeLength to targetEdgeLength
      edgeCenter = lerp(verts[v0],verts[v1],.5)
      newVerts[v0] = vpv(newVerts[v0], lerp(edgeCenter, verts[v0], targetEdgeLength/currentEdgeLength))
      newVerts[v1] = vpv(newVerts[v1], lerp(edgeCenter, verts[v1], targetEdgeLength/currentEdgeLength))
      nContributions[v0] += 1
      nContributions[v1] += 1
    for i in range(len(newVerts)):
      if nContributions[i] == 0:
        newVerts[i] = verts[i]
      else:
        newVerts[i] = vxs(newVerts[i], 1./nContributions[i])
    verts = newVerts

    # Project (exactly) to the average-radius sphere.
    radii = [length(vert) for vert in verts]
    avgRadius = sum(radii) / len(radii)
    verts = [vxs(vert, avgRadius/radius) for vert,radius in zip(verts,radii)]
  t1 = time.time()

  if verboseLevel >= 1: sys.stderr.write("\n")
  if verboseLevel >= 1: sys.stderr.write("      did %r iterations in %.6f secs.\n"%(nIters, t1-t0))
  if verboseLevel >= 1: sys.stderr.write("      final verts = %r\n"%(verts,))
  if verboseLevel >= 1: sys.stderr.write("      final spring lengths with '==' constraints: %r\n"%([dist(verts[v0],verts[v1]) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation=="=="]))
  if verboseLevel >= 1: sys.stderr.write("      final spring lengths with '>=' constraints: %r\n"%([dist(verts[v0],verts[v1]) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation==">="]))
  if verboseLevel >= 1: sys.stderr.write("      final spring lengths with '<=' constraints: %r\n"%([dist(verts[v0],verts[v1]) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation=="<="]))
  if verboseLevel >= 1: sys.stderr.write("      final radii = %r\n"%([length(vert) for vert in verts]))
  sys.stderr.write("    out relaxOnAnySizeSphere\n")
  return verts


# Show three projections of a graph scaled to the unit sphere.
def printSVG(verts, edges):
  sys.stderr.write("verts = %r\n"%(verts,))
  sys.stderr.write("edges = %r\n"%(edges,))
  print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>")
  print("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"768pt\" height=\"128pt\" viewBox=\"0 0 768 128\" >")

  maxRadius = max(length(vert) for vert in verts)
  scale = 50/maxRadius  # can't put this as part of the xform, or it will thicken the lines
  pointRadius = .02

  # show three projections of it.
  translatex = 64  # for starters
  for xy in [[0,1], [0,2], [1,2]]:  # xy, xz, yz.  no particular reason except they look ok for frucht.
    print(" <g style=\"fill:none;stroke:#000000\" transform=\"translate(%r,64)\" >"%(translatex,))
    for vert in verts:
      x = vert[xy[0]]
      y = vert[xy[1]]
      print("  <circle cx=\"%r\" cy=\"%r\" r=\"%r\" />"%(x*scale,-y*scale,pointRadius*scale))
    for v0,v1 in edges:
      x0 = verts[v0][xy[0]]
      y0 = verts[v0][xy[1]]
      x1 = verts[v1][xy[0]]
      y1 = verts[v1][xy[1]]
      print("  <path d=\"M %r,%r L %r,%r\" />"%(x0*scale,-y0*scale,x1*scale,-y1*scale))
    print(" </g>")
    translatex += 128
    
  print("</svg>")


if len(sys.argv) != 2:
  exit("Usage: relax.py <nIters> < input.txt > output.svg")
nIters = int(sys.argv[1])

verts = []
edgeLengthConstraints = []

if sys.stdin.isatty():
  sys.stderr.write("Enter 3d vertex coords (3 numbers per line), and distance constraints:")
  sys.stderr.flush()
for line in sys.stdin:
  # strip optional comment (beginning with '#') off the end, and whitespace
  line = re.sub('#.*', '', line).strip()
  if line == '': continue
  #sys.stderr.write("line = %r\n"%(line))
  tokens = line.split()
  if len(tokens) == 3:
    # vertex coords: x y z
    verts.append([float(token) for token in tokens])
  elif len(tokens) == 4:
    # it's an edge length constraint, e.g. "3 6 >= 1.23"
    v0 = int(tokens[0])
    v1 = int(tokens[1])
    relation = tokens[2]
    targetEdgeLength = float(tokens[3])
    if relation not in ['==','<=','>=']:
      exit("Bad line %r"%(line,))
    edgeLengthConstraints.append([v0,v1,relation,targetEdgeLength])
  else:
    exit("bad line %r" % (line,))


verts = relaxOnAnySizeSphere(verts, edgeLengthConstraints, nIters)

# Assume the actual edges are the springs with target exactly 1.
edges = [(v0,v1) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation=="==" and targetEdgeLength==1]
printSVG(verts, edges)

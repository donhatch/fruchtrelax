#!/usr/bin/python3

import math
import re
import sys

if sys.stdin.isatty():
  print("Enter 3d vertex coords, separated by spaces and/or newlines:")

verboseLevel = 2  # hard coded.  0: just output, 1: reasonable interesting progress messages, 2: debug


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
  verboseLevel = 2  # manually set this to something higher to debug
  if verboseLevel >= 1: print("    in relaxOnAnySizeSphere")
  if verboseLevel >= 1: print("      initial verts = %r"%(verts,))
  if verboseLevel >= 1: print("      edgeLengthConstraints = %r"%(edgeLengthConstraints,))
  if verboseLevel >= 1: print("      nIters = %r"%(nIters,))

  for iIter in range(nIters):
    if verboseLevel >= 2: print("      =============================")
    if verboseLevel >= 2: print("      iIter=%d/%d"%(iIter,nIters))
    if verboseLevel >= 2: print("          before:")
    if verboseLevel >= 3: print("              verts = %r"%(verts,))
    if verboseLevel >= 2: print("              edgeLengths = %r"%([dist(verts[v0],verts[v1]) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation=="=="]))
    if verboseLevel >= 2: print("              radii = %r"%([length(vert) for vert in verts]))

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
      # this edge wants to change (or keep) its length from currentEdgeLength to targetEdgeLength
      if verboseLevel >= 3: print("                edge (%r,%r) wants to change length from %r to %r"%(v0,v1,currentEdgeLength,targetEdgeLength))
      edgeCenter = lerp(verts[v0],verts[v1],.5)
      newVerts[v0] = vpv(newVerts[v0], lerp(edgeCenter, verts[v0], targetEdgeLength/currentEdgeLength))
      newVerts[v1] = vpv(newVerts[v1], lerp(edgeCenter, verts[v1], targetEdgeLength/currentEdgeLength))
      nContributions[v0] += 1
      nContributions[v1] += 1
    if verboseLevel >= 3: print("          newVerts (before divide) = %r"%(newVerts,))
    if verboseLevel >= 3: print("          nContributions = %r"%(nContributions,))
    for i in range(len(newVerts)):
      if nContributions[i] == 0:
        newVerts[i] = verts[i]
      else:
        newVerts[i] = vxs(newVerts[i], 1./nContributions[i])
    verts = newVerts
    if verboseLevel >= 3: print("          verts = %r"%(verts,))

    if verboseLevel >= 3: print("          after edge relax, before sphere project:")
    if verboseLevel >= 3: print("              verts = %r"%(verts,))
    if verboseLevel >= 3: print("              edgeLengths = %r"%([dist(verts[v0],verts[v1]) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation=="=="]))
    if verboseLevel >= 3: print("              radii = %r"%([length(vert) for vert in verts]))

    # Project (exactly) to average-radius sphere
    radii = [length(vert) for vert in verts]
    avgRadius = sum(radii) / len(radii)
    print("              avgRadius = %r"%(avgRadius,))
    verts = [vxs(vert, avgRadius/radius) for vert,radius in zip(verts,radii)]
    if verboseLevel >= 3: print("          after:")
    if verboseLevel >= 3: print("              verts = %r"%(verts,))
    if verboseLevel >= 3: print("              edgeLengths = %r"%([dist(verts[v0],verts[v1]) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation=="=="]))
    if verboseLevel >= 3: print("              radii = %r"%([length(vert) for vert in verts]))
    print("      =============================")
  print("    out relaxOnAnySizeSphere")
  return verts


def printSVG(verts, edges):
  sys.stderr.write("verts = %r\n"%(verts,))
  sys.stderr.write("edges = %r\n"%(edges,))
  print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>")
  #print("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"256pt\" height=\"256pt\" viewBox=\"0 0 256 256\" transform=\"translate(128,128)\">")
  print("<svg xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:cc=\"http://creativecommons.org/ns#\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"768pt\" height=\"256pt\" viewBox=\"0 0 768 256\" >")

  maxRadius = max(length(vert) for vert in verts)
  scale = 50/maxRadius  # can't put this as part of the xform, or it will thicken the lines
  pointRadius = .02

  print(" <g style=\"fill:none;stroke:#000000\" transform=\"translate(128,128)\" >")
  for x,y,z in verts:
    print("  <circle cx=\"%r\" cy=\"%r\" r=\"%r\" />"%(x*scale,-y*scale,pointRadius*scale))
  for v0,v1 in edges:
    x0,y0,z0 = verts[v0]
    x1,y1,z1 = verts[v1]
    print("  <path d=\"M %r,%r L %r,%r\" />"%(x0*scale,-y0*scale,x1*scale,-y1*scale))
  print(" </g>")

  print(" <g style=\"fill:none;stroke:#000000\" transform=\"translate(256,128)\" >")
  for x,y,z in verts:
    print("  <circle cx=\"%r\" cy=\"%r\" r=\"%r\" />"%(y*scale,-z*scale,pointRadius*scale))
  for v0,v1 in edges:
    x0,y0,z0 = verts[v0]
    x1,y1,z1 = verts[v1]
    print("  <path d=\"M %r,%r L %r,%r\" />"%(y0*scale,-z0*scale,y1*scale,-z1*scale))
  print(" </g>")

  print(" <g style=\"fill:none;stroke:#000000\" transform=\"translate(384,128)\" >")
  for x,y,z in verts:
    print("  <circle cx=\"%r\" cy=\"%r\" r=\"%r\" />"%(x*scale,-z*scale,pointRadius*scale))
  for v0,v1 in edges:
    x0,y0,z0 = verts[v0]
    x1,y1,z1 = verts[v1]
    print("  <path d=\"M %r,%r L %r,%r\" />"%(x0*scale,-z0*scale,x1*scale,-z1*scale))
  print(" </g>")

  print("</svg>")


#
# .2352 .234 .2345
# .3456435 .3456435 .3567
# 0 1 == .33333
# 0 2 <= .33333
# 0 3 >= .33333
if True:
  if len(sys.argv) != 2:
    exit("Usage: relax.py <nIters>")
  nIters = int(sys.argv[1])

  verts = []
  edgeLengthConstraints = []
  for line in sys.stdin:
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

  edges = [(v0,v1) for v0,v1,relation,targetEdgeLength in edgeLengthConstraints if relation=="==" and targetEdgeLength==1]
  printSVG(verts, edges)

import numpy as np
import math as m
import PlotLoop

IMAX = 11 # grid poins are i = 0, ..., IMAX - 1
L = 1.0
DOF = 3
C = 1.0
DT = 0.005
ITMAX = 500 
DEBUG = True
AXISRANGE = 0.0, 1.0, -2.0, 2.0

# i    i+1   i+2   i+3  <-- face
# +-----+-----+-----+
#    i    i+1   i+2     <-- cell

# x_i   x_i+1  <-- x
# -+------+-
# -1      1    <-- xi

def GaussianQuadrature(n):
	# Gaussian quadrature points and weights for integral [-1, 1]

	if n == 1:
		return np.array([0.0]), np.array([2.0])
	elif n == 2:
		return np.array([-m.sqrt(3) / 3.0, m.sqrt(3) / 3.0]), np.array([1.0, 1.0])
	elif n == 3:
		return np.array([-m.sqrt(0.6), 0.0, m.sqrt(0.6)]), np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
	else:
		raise ValueError

def Flux(u):

	return C * u

def NumericalFlux(ul, ur):

	return C * 0.5 * (ul + ur) - 0.5 * abs(C) * (ur - ul)

class LogicalTaylorBasis(object):

	def __init__(self):
		self.Order = 2
		pass

	def Evaluate(self, xi):
		return np.array([1.0, xi, 0.5 * xi * xi])

	def EvaluateDerivative(self, xi):
		return np.array([0.0, 1.0, xi])

class Cell1D(object):
	""" xi = -1.0 -> 1.0 """

	def __init__(self, icell):

		self.dx = L / (IMAX - 1)
		self.ID = icell
		self.basis = LogicalTaylorBasis()

		n = self.basis.Order
		o1 = n * n # polynomial order for mass matrix
		o2 = n * (n - 1) # polynomial order of flux-volume integral
		self.ngauss1 = int((o1 + 1) / 2.0 + 1.0)
		self.ngauss2 = int((o2 + 1) / 2.0 + 1.0)

	def Jacobian(self, xi):

		return 0.5 * self.dx

	def EvaluateCellU(self, xi, U):

		phi = self.basis.Evaluate(xi)
		u = np.dot(U[self.ID, :], phi)
		return u, phi

	def IntegrateFluxVolume(self, R, U):

		xis, ws = GaussianQuadrature(self.ngauss2)

		for k in range(len(xis)):
			xi, w = xis[k], ws[k]
			J = self.Jacobian(xi)
			u, phi = self.EvaluateCellU(xi, U)
			f = Flux(u)
			dphidxi = self.basis.EvaluateDerivative(xi)
			dphidx = dphidxi / J # FIXME: verify this. probably only valid for 1-d problems
			R[self.ID, :] += w * f * dphidx * J
		#print "FluxVolume %02d: " % self.ID, R[self.ID, :]

	def MassMatrix(self):

		xis, ws = GaussianQuadrature(self.ngauss1)
		M = np.zeros((6))
		for k in range(len(xis)):
			xi, w = xis[k], ws[k]
			J = self.Jacobian(xi)
			phi = self.basis.Evaluate(xi)
			M[0] += w * phi[0] * phi[0] * J # M11
			M[1] += w * phi[0] * phi[1] * J # M12
			M[2] += w * phi[1] * phi[1] * J # M22
			M[3] += w * phi[0] * phi[2] * J # M13
			M[4] += w * phi[1] * phi[2] * J # M23
			M[5] += w * phi[2] * phi[2] * J # M33
		return M

class Face1D(object):

	def __init__(self, cellL, cellR):

		self.L = cellL
		self.R = cellR

	def IntegrateFluxSurface(self, R, U):

		labelL, labelR = "No", "No"
		if self.L != None:
			uL, phiL = self.L.EvaluateCellU(1.0, U)
			labelL = "%02d" % self.L.ID
		if self.R != None:
			uR, phiR = self.R.EvaluateCellU(-1.0, U)
			labelR = "%02d" % self.R.ID
		label = "%s-%s" % (labelL, labelR)

		if self.L == None:
			uL = uR
			phiL = None
		if self.R == None:
			uR = uL 
			phiR = None

		f = NumericalFlux(uL, uR)
		print "NumFlux ", label , " uL, uR = ", uL, uR, ", f = ", f, "phiL, phiR = ", phiL, phiR

		if self.L != None:
			print "Add to Cell ", self.L.ID, -f * phiL
			R[self.L.ID, :] -= f * phiL
		if self.R != None:
			print "Add to Cell ", self.R.ID, f * phiR
			R[self.R.ID, :] += f * phiR

def EvaluateResidual(R, U, cells, faces):

	if DEBUG:
		print "R: init"
		print R
	for cell in cells:
		cell.IntegrateFluxVolume(R, U)
	if DEBUG:
		print "R: after volume integral"
		print R

	for face in faces:
		face.IntegrateFluxSurface(R, U)
	if DEBUG:
		print "R: after surface integral"
		print R

def EvaluateInverseMassMatrix(Minv, cells):

	for cell in cells:
		m = cell.MassMatrix()
		m = np.matrix(FromPackedToFull(m))
		Minv[cell.ID, :] = FromFullToPacked(m.I)

def MinvR(Minv, R):
	""" applies Minv to R """

	for i in range(Minv.shape[0]):
		minv = FromPackedToFull(Minv[i, :])
		print "minv: %02d" % i, minv
		R[i, :] = np.dot(minv, R[i, :])

class ExplicitIntegrator(object):

	def __init__(self, cells, faces, U):

		N = len(cells)
		Minv = np.zeros((N, DOF * (DOF + 1) / 2))
		R = np.zeros((N, DOF))
		dt = DT

		EvaluateInverseMassMatrix(Minv, cells)
		print "Minv"
		print Minv
		self.Minv = Minv
		self.U = U
		self.R = R
		self.dt = dt

		self.Cells = cells
		self.Faces = faces

		self.Iteration = 0

	def Step(self):

		self.R[:, :] = 0.0
		EvaluateResidual(self.R, self.U, self.Cells, self.Faces)
		MinvR(self.Minv, self.R)
		print "DT = ", self.dt

		self.U[:, :] += self.dt * self.R

		self.Iteration += 1

def DumpCellU(f, cells, U):

	dx = L / float(IMAX - 1)
	for cell in cells:
		for xi in np.linspace(-1.0, 1.0, 10):
			u, phi = cell.EvaluateCellU(xi, U)
			x = float(cell.ID) * dx + (xi + 1.) / 2.0 * dx
			print >>f, x, u

def ReconstructCellUs(cells, U):

	NR = 10 # # of points used within each cell to represent reconstructed solution

	numCells = len(cells)

	X = np.zeros((numCells, NR + 1)) # +1 for inserting np.nan, needed to show discontinuity at inter-cell boundaries
	UU = np.zeros((numCells, NR + 1)) # +1 for inserting np.nan, needed to show discontinuity at inter-cell boundaries

	dx = L / float(IMAX - 1)
	for i, cell in enumerate(cells):
		X[i, 0] = np.nan
		UU[i, 0] = np.nan
		for j, xi in enumerate(np.linspace(-1.0, 1.0, NR)):
			u, phi = cell.EvaluateCellU(xi, U)
			x = float(cell.ID) * dx + (xi + 1.) / 2.0 * dx
			X[i, j + 1] = x
			UU[i, j + 1] = u

	X = X.reshape(numCells * (NR + 1))
	UU = UU.reshape(numCells * (NR + 1))

	return X, UU

def SineWave(x):
	a = 1.0
	dx = L / (IMAX - 1)
	J = 0.5 * dx
	return (
		m.sin(a * 2.0 * m.pi * x),
		a * 2.0 * m.pi * m.cos(a * 2.0 * m.pi * x) * J,
		-a * a * 4.0 * m.pi * m.pi * m.sin(a * 2.0 * m.pi * x) * J * J
		)

def SineWave2(x):
	a = 1.0
	b = 4.0
	A = 1.0
	B = 0.2
	dx = L / (IMAX - 1)
	J = 0.5 * dx
	a2pi = a * 2.0 * m.pi
	b2pi = b * 2.0 * m.pi
	return (
		A * m.sin(a2pi * x) + B * m.sin(b2pi * x),
		(A * a2pi * m.cos(a2pi * x) + B * b2pi * m.cos(b2pi * x)) * J,
		(-A * a2pi * a2pi * m.sin(a2pi * x) - B * b2pi * b2pi * m.sin(b2pi * x)) * J * J
		)

def Constant(x):
	return 1.0, 0.0, 0.0

def FromPackedToFull(m):

	return np.array([
		[m[0], m[1], m[3]],
		[m[1], m[2], m[4]],
		[m[3], m[4], m[5]]
		])

def FromFullToPacked(m):

	return np.array([m[0, 0], m[0, 1], m[1, 1], m[0, 2], m[1, 2], m[2, 2]])

def Plot(cells, U):

	XX, UU = ReconstructCellUs(cells, U)
	PlotLoop.Clear()
	PlotLoop.GetPlot().plot(XX, UU)
	PlotLoop.GetPlot().axis(AXISRANGE)
	PlotLoop.Draw()

def Init(U):

	#f = Constant
	#f = SineWave
	f = SineWave2
	for i in range(0, IMAX - 1):
		x1 = float(i) / float(IMAX - 1)
		x2 = float(i + 1) / float(IMAX - 1)
		x = 0.5 * (x1 + x2)
		u, up, upp = f(x)
		U[i, :] = np.array([u, up, upp])

def Step(integrator, cells, U):

	integrator.Step()
	Plot(cells, U)

	if integrator.Iteration >= ITMAX:
		return False

	return True

def Main():

	PlotLoop.Init()

	dx = 1.0 / float(IMAX - 1)
	cfl = C * DT / dx
	print "CFL = ", cfl

	U = np.zeros((IMAX - 1, DOF))
	Init(U)

	cells = [ Cell1D(i) for i in range(IMAX - 1) ]

	faces = []
	faces.append(Face1D(cells[-1], cells[0])) # periodic bounday at x = 0.0
	for ic in range(1, len(cells)):
		faces.append(Face1D(cells[ic - 1], cells[ic]))
	# periodic boundary at x = L is already taken care of by the face at x = 0.0

	integrator = ExplicitIntegrator(cells, faces, U)

	Plot(cells, U)

	PlotLoop.RegisterLoopCallback(lambda : Step(integrator, cells, U))

	"""
	for i in range(ITMAX):
		integrator.Step()
		for ic in range(IMAX - 1):
			print "Cell %d" % ic, U[ic, :]
		#DumpU(open("u.%03d.dat" % i, "w"), U)
		DumpCellU(open("u.%03d.dat" % i, "w"), cells, U)
	"""

	PlotLoop.Start()

if __name__ == "__main__":
	Main()


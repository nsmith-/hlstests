import math
from operator import mul

max_iter = 24

def format(a, w=4, ind="  "):
    for i in range(len(a)/w):
        print "%s%s," % (ind, ",".join(str(x) for x in a[i*w:(i+1)*w]), )


print "#ifndef LORENTZ_TABLES_H"
print "#define LORENTZ_TABLES_H"

ap_lorentz_cordic_angles = [math.atan(pow(2, -i-1)) for i in range(max_iter)]
print "const double ap_lorentz_cordic_angles[%d] = {" % max_iter
format(ap_lorentz_cordic_angles, 6)
print "};"

iter_scales = [pow(1+pow(2, -2*(i+1)), -.5) for i in range(max_iter)]
ap_lorentz_cordic_scales = [reduce(mul, iter_scales[:n+1]) for n in range(max_iter)]
print "const double ap_lorentz_cordic_scales[%d] = {" % max_iter
format(ap_lorentz_cordic_scales, 6)
print "};"

ap_lorentz_cordic_angles_hyp = [math.atanh(pow(2, -i-2)) for i in range(max_iter)]
print "const double ap_lorentz_cordic_angles_hyp[%d] = {" % max_iter
format(ap_lorentz_cordic_angles_hyp, 6)
print "};"

iter_scales_hyp = [pow(1-pow(2, -2*(i+2)), -.5) for i in range(max_iter)]
lutsize = 16

ap_lorentz_cordic_xexact_hyp = [math.cosh(eta*0.5) for eta in range(lutsize)]
print "const double ap_lorentz_cordic_xexact_hyp[%d] = {" % lutsize
format(ap_lorentz_cordic_xexact_hyp, 4)
print "};"

ap_lorentz_cordic_yexact_hyp = [math.sinh(eta*0.5) for eta in range(lutsize)]
print "const double ap_lorentz_cordic_yexact_hyp[%d] = {" % lutsize
format(ap_lorentz_cordic_yexact_hyp, 4)
print "};"

ap_lorentz_cordic_xcoarse_hyp = [ [math.cosh(eta*0.5)*reduce(mul, iter_scales_hyp[:n+1]) for eta in range(lutsize)] for n in range(max_iter)]
print "const double ap_lorentz_cordic_xcoarse_hyp[%d][%d] = {" % (max_iter, lutsize)
for i in range(max_iter):
    print "  {"
    format(ap_lorentz_cordic_xcoarse_hyp[i], 4, "    ")
    print "  }," if i<max_iter-1 else "  }"
print "};"

ap_lorentz_cordic_ycoarse_hyp = [ [math.sinh(eta*0.5)*reduce(mul, iter_scales_hyp[:n+1]) for eta in range(lutsize)] for n in range(max_iter)]
print "const double ap_lorentz_cordic_ycoarse_hyp[%d][%d] = {" % (max_iter, lutsize)
for i in range(max_iter):
    print "  {"
    format(ap_lorentz_cordic_ycoarse_hyp[i], 4, "    ")
    print "  }," if i<max_iter-1 else "  }"
print "};"

print "#endif // LORENTZ_TABLES_H"

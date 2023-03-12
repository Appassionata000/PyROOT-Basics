import math
import ROOT as r

def double_sided_crystal_ball(x, p):
    # parameters: p[0]: constants, p[1]: mean, p[2]: std 
    #             p[3]: a1(low),  p[4]: n1(low),  p[5]: a2(high)  p[6]: n2(high)
    s = (x[0] - p[1]) / p[2]
    a1, n1, a2, n2 = p[3], p[4], p[5], p[6]
    
    result = 1
        
    if (n1 / abs(a1) - abs(a1) - s) > 0 and (n2 / abs(a2) - abs(a2) + s) > 0:

        A1 = math.pow((n1 / abs(a1)), n1) * math.exp(-a1**2 / 2)
        A2 = math.pow((n2 / abs(a2)), n2) * math.exp(-a2**2 / 2)
        B1 = math.pow((n1 / abs(a1) - abs(a1) - s), -n1)
        B2 = math.pow((n2 / abs(a2) - abs(a2) + s), -n2)

        if s < -a1:
            result = p[0] * A1 * B1
        elif s > a2:
            result = p[0] * A2 * B2
        else:
            result = p[0] * math.exp(-s**2 / 2)

    return result

def cubic_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):
    
        r.TF1.RejectPoint(True)
    
    return p[0] + p[1] * x[0] + p[2] * x[0]**2 + p[3] * x[0]**3 

def quartic_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):
    
        r.TF1.RejectPoint(True)
        
    return p[0] + p[1] * x[0] + p[2] * x[0]**2 + p[3] * x[0]**3 + p[4] * x[0]**4

def quintic_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):
    
        r.TF1.RejectPoint(True)
        
    return p[0] + p[1] * x[0] + p[2] * x[0]**2 + p[3] * x[0]**3 + p[4] * x[0]**4 + p[5] * x[0]**5

def sextic_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):
    
        r.TF1.RejectPoint(True)
        
    return p[0] + p[1] * x[0] + p[2] * x[0]**2 + p[3] * x[0]**3 + p[4] * x[0]**4 + p[5] * x[0]**5 + p[6] * x[0]**6

def expPoly1_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):

        r.TF1.RejectPoint(True)
    
    return math.exp(p[0] * x[0])

def expPoly2_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):

        r.TF1.RejectPoint(True)
    
    return p[0] * math.exp(-p[1] * x[0] - p[2] * x[0]**2)

def expPoly3_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):

        r.TF1.RejectPoint(True)
    
    return p[0] * math.exp(-p[1] * x[0] - p[2] * x[0]**2 - p[3] * x[0]**3)

def expPoly4_bkg(x, p):
    
    if (x[0] >= 120e3 and x[0] <= 130e3):
        
        r.TF1.RejectPoint(True)
        
    return p[0] * math.exp(-p[1] * x[0] - p[2] * x[0]**2 - p[3] * x[0]**3 - p[4] * x[0]**4)

def gaus_1(x, p):
    
    return p[0] * math.exp(-0.5 * ((x[0] - p[1]) / p[2])**2)
    
def gaus_3(x, p):
    
    gaus1 = p[0] * math.exp(-0.5 * ((x[0] - p[1]) / p[2])**2)
    gaus2 = p[3] * math.exp(-0.5 * ((x[0] - p[4]) / p[5])**2)
    gaus3 = p[6] * math.exp(-0.5 * ((x[0] - p[7]) / p[8])**2)
   
    return gaus1 + gaus2 + gaus3


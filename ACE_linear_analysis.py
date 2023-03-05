from gurobipy import *

def Linear_Constr_branch(Var, dummy):
    Constr = []
    Constr.append(' + '.join(Var) + ' - 2 ' + dummy + ' = 0')
    return Constr

#Var表示乘数，pro表示乘积，都是二进制
def Constr_mul(Var, pro):
    Constr = []
    len_Var = len(Var)
    con = ' + '.join(Var)
    Constr.append(con + ' - ' + str(len_Var) + ' ' + pro + ' >= 0')
    Constr.append(con + ' - ' + str(len_Var) + ' ' + pro + ' <= ' + str(len_Var - 1))
    return Constr

#pro = (1-var1)*var2*var3*(1-var4), 都是二进制
def Constr_mul4(Var, pro):
    Constr = []
    con = Var[1] + ' + ' + Var[2] + ' - ' + Var[0] + ' - ' + Var[3]
    Constr.append(con + ' - 4 ' + pro + ' >= -2')
    Constr.append(con + ' - 2 ' + pro + ' <= 1')
    return Constr

#其中Im, L1, L2, Om, Beta是n维2进制变量，beta是一个2进制变量，Dummy是n维整数变量，dummy1,dummy2是一个整数变量
#L1[i]表示第i个与门左边的输入,L2[i]表示第i个与门右边的输入
def Linear_Constr_g(Im, L1, L2, Om, Beta, beta, Dummy, dummy1, dummy2):
    Constr = []
    obj1 = []
    obj2 = []
    obj3 = []

    n = len(Im)
    for i in range(n):
        #(4)
        Constr = Constr + Linear_Constr_branch([Im[i], L1[i], L2[(i-1)%n]], Dummy[i])

        #(6)
        Constr = Constr + [L1[i] + ' - ' + Om[i] + ' <= 0']
        Constr = Constr + [L2[i] + ' - ' + Om[i] + ' <= 0']

        #(7)
        Constr = Constr + Constr_mul4([Beta[(i-1)%n], Om[i], Om[(i+1)%n], beta], Beta[i])
        Constr = Constr + [L1[i] + ' - ' + L2[(i + 1) % n] + ' + ' + Beta[i] + ' <= 1']
        Constr = Constr + [L2[(i + 1) % n] + ' - ' + L1[i] + ' + ' + Beta[i] + ' <= 1']

    #(5)
    Constr = Constr + Constr_mul(Om, beta)

    #(8)
    N = int(n / 2)
    set1 = [Im[2*i] for i in range(N)]
    set2 = [Im[2*i+1] for i in range(N)]

    Constr = Constr + [' + '.join(set1) + ' + ' + dummy1 + ' + [ ' + beta + ' * ' + dummy1 + ' ] = 0']
    Constr = Constr + [' + '.join(set2) + ' + ' + dummy2 + ' + [ ' + beta + ' * ' + dummy2 + ' ] = 0']

    obj1 = obj1 + Om
    obj2 = obj2 + Beta
    obj3 = obj3 + [beta]

    return [Constr, obj1, obj2, obj3]

def Linear_Constr_ACE(Step, PI):
    obj1 = []
    obj2 = []
    obj3 = []
    Constr = []
    bin_var = []
    general_var = []
    inv_PI = [0 for i in range(5)]
    for i in range(5):
        inv_PI[PI[i]] = i

    A = [['A_s' + str(j) + '_' + str(i) for i in range(64)] for j in range(Step+1)]
    B = [['B_s' + str(j) + '_' + str(i) for i in range(64)] for j in range(Step+1)]
    C = [['C_s' + str(j) + '_' + str(i) for i in range(64)] for j in range(Step+1)]
    D = [['D_s' + str(j) + '_' + str(i) for i in range(64)] for j in range(Step+1)]
    E = [['E_s' + str(j) + '_' + str(i) for i in range(64)] for j in range(Step+1)]
    S = [A,B,C,D,E]
    for i in range(Step+1):
        for j in range(64):
            bin_var.append(A[i][j])
            bin_var.append(B[i][j])
            bin_var.append(C[i][j])
            bin_var.append(D[i][j])
            bin_var.append(E[i][j])

    #额外的限制条件
    '''
    #1
    SR1 = []
    SR2 = []
    for i in range(32, 64):
        SR1.append(A[0][i])
        SR1.append(C[0][i])
        SR2.append(A[Step][i])
        SR2.append(C[Step][i])
    Constr.append(' + '.join(SR1) + ' >= 1')
    Constr.append(' + '.join(SR2) + ' >= 1')

    SC1 = []
    SC2 = []
    for i in range(32):
        SC1.append(A[0][i])
        SC1.append(C[0][i])
        SC2.append(A[Step][i])
        SC2.append(C[Step][i])
    for i in range(64):
        SC1.append(B[0][i])
        SC1.append(D[0][i])
        SC1.append(E[0][i])
        SC2.append(B[Step][i])
        SC2.append(D[Step][i])
        SC2.append(E[Step][i])
    Constr.append(' + '.join(SC1) + ' = 0')
    Constr.append(' + '.join(SC2) + ' = 0')
    '''

    #2
    IN = []
    OUT = []
    for i in range(5):
        for j in range(64):
            IN.append(S[i][0][j])
            OUT.append(S[i][Step][j])

    Constr.append(' + '.join(IN) + ' >= 1')
    Constr.append(' + '.join(OUT) + ' >= 1')


    gimA = [[['gimA_s' + str(i) + '_r' + str(j) + '_' + str(k) for k in range(32)] for j in range(8)] for i in range(Step)]
    gomA = [[['gomA_s' + str(i) + '_r' + str(j) + '_' + str(k) for k in range(32)] for j in range(8)] for i in range(Step)]
    gimC = [[['gimC_s' + str(i) + '_r' + str(j) + '_' + str(k) for k in range(32)] for j in range(8)] for i in range(Step)]
    gomC = [[['gomC_s' + str(i) + '_r' + str(j) + '_' + str(k) for k in range(32)] for j in range(8)] for i in range(Step)]
    gimE = [[['gimE_s' + str(i) + '_r' + str(j) + '_' + str(k) for k in range(32)] for j in range(8)] for i in range(Step)]
    gomE = [[['gomE_s' + str(i) + '_r' + str(j) + '_' + str(k) for k in range(32)] for j in range(8)] for i in range(Step)]
    gim = [gimA, gimC, gimE]
    gom = [gomA, gomC, gomE]
    x8A = [['x8A_s' + str(i) + '_' + str(j) for j in range(32)] for i in range(Step)]
    x8C = [['x8C_s' + str(i) + '_' + str(j) for j in range(32)] for i in range(Step)]
    x8E = [['x8E_s' + str(i) + '_' + str(j) for j in range(32)] for i in range(Step)]
    x8 = [x8A, x8C, x8E]

    for i in range(3):
        for j in range(Step):
            for k in range(8):
                for l in range(32):
                    bin_var.append(gim[i][j][k][l])
                    bin_var.append(gom[i][j][k][l])

    for i in range(3):
        for j in range(Step):
            for k in range(32):
                bin_var.append(x8[i][j][k])

    #ACE线性刻画
    for i in range(Step):

        Sout = [S[inv_PI[0]][i+1], S[inv_PI[1]][i+1], S[inv_PI[2]][i+1], S[inv_PI[3]][i+1], S[inv_PI[4]][i+1]]

        #对分支的刻画：
        dummyA = ['dummyA_s' + str(i) + '_' + str(j) for j in range(64)]
        dummyC = ['dummyC_s' + str(i) + '_' + str(j) for j in range(64)]
        dummyE = ['dummyE_s' + str(i) + '_' + str(j) for j in range(64)]
        for j in range(64):
            general_var.append(dummyA[j])
            general_var.append(dummyC[j])
            general_var.append(dummyE[j])
        branchA = x8A[i] + gomA[i][7]
        branchC = x8C[i] + gomC[i][7]
        branchE = x8E[i] + gomE[i][7]

        for j in range(64):
            Constr = Constr + Linear_Constr_branch([branchA[j], Sout[4][j], Sout[0][j]], dummyA[j])
            Constr = Constr + Linear_Constr_branch([branchC[j], Sout[1][j], Sout[2][j]], dummyC[j])
            Constr = Constr + Linear_Constr_branch([branchE[j], Sout[3][j], Sout[4][j]], dummyE[j])

        #异或刻画
        for j in range(64):
            Constr = Constr + [S[1][i][j] + ' - ' + Sout[1][j] + ' = 0']
            Constr = Constr + [S[3][i][j] + ' - ' + Sout[3][j] + ' = 0']

        #SB-64刻画
        for j in range(3):
            for k in range(8):
                #刻画环型与门组合
                Gim = [None for l in range(32)]
                Gom = [None for l in range(32)]
                for l in range(32):
                    Gim[l] = gim[j][i][k][(-l*5)%32]
                    Gom[l] = gom[j][i][k][(-l*5)%32]

                L1 = ['L1' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                L2 = ['L2' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                Beta = ['Beta' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                beta = 'beta' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k)
                Dummyg = ['Dummyg' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                dummyg1 = 'dummyg1' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k)
                dummyg2 = 'dummyg2' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k)

                for l in range(32):
                    bin_var.append(L1[l])
                    bin_var.append(L2[l])
                    bin_var.append(Beta[l])
                    general_var.append(Dummyg[l])
                bin_var.append(beta)
                general_var.append(dummyg1)
                general_var.append(dummyg2)

                [temp_Constr, temp_obj1, temp_obj2, temp_obj3] = Linear_Constr_g(Gim, L1, L2, Gom, Beta, beta, Dummyg, dummyg1, dummyg2)
                Constr = Constr + temp_Constr
                obj1 = obj1 + temp_obj1
                obj2 = obj2 + temp_obj2
                obj3 = obj3 + temp_obj3

                #刻画SB-64轮函数中的分支
                dummysr = ['dummysr' + chr(j*2 + 65) + '_s' + str(i) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                general_var = general_var + dummysr

                branch2 = gim[j][i][k]
                branch3 = [None for l in range(32)]
                for l in range(32):
                    branch3[l] = gom[j][i][k][(l+1)%32]

                if k == 0:
                    branch1 = S[j*2][i][32:64]
                else:
                    branch1 = gom[j][i][k-1]
                if k == 7:
                    branch4 = x8[j][i]
                else:
                    branch4 = gom[j][i][k+1]

                for l in range(32):
                    Constr = Constr + Linear_Constr_branch([branch1[l], branch2[l], branch3[l], branch4[l]], dummysr[l])

                #刻画SB-64轮函数中的相等
                for l in range(32):
                    Constr = Constr + [S[j*2][i][l] + ' - ' + gom[j][i][0][l] + ' = 0']

    return [obj1, obj2, obj3, Constr, bin_var, general_var]

def ACELinearmodel_Generate(Step, PI):
    [obj1, obj2, obj3, Constr, bin_var, general_var] = Linear_Constr_ACE(Step, PI)

    filename = 'ACELinearmodel_S' + str(Step) + '.lp'
    file = open(filename, 'w')

    #objective
    Obj = ' + '.join(obj1)
    for o in obj2:
        Obj = Obj + ' - ' + o
    for o in obj3:
        Obj = Obj + ' - 17 ' + o
    file.write('min\n' + Obj + '\n')

    file.write('Subject to\n')
    for c in Constr:
        file.write(c + '\n')

    file.write('Binary\n')
    for bv in bin_var:
        file.write(bv + '\n')

    file.write('General\n')
    for gv in general_var:
        file.write(gv + '\n')

    file.close()

def mycallback(model, where):
    if where == GRB.Callback.MIPSOL:
        CB = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
        if CB == 212:
            m.terminate()

if __name__ == '__main__':
    PI = [3,2,0,4,1]
    Step = 8
    ACELinearmodel_Generate(Step, PI)
    modelname = 'ACELinearmodel_S' + str(Step)
    m = read(modelname)
    m.setParam(GRB.Param.MIPFocus, 2)
    m.optimize(mycallback)
    vars = m.getVars()
    File = open('Result_ACELinear_S' + str(Step) + '.txt', 'w')
    for v in vars:
        if (v.varName[0:2] in ['A_', 'B_', 'C_', 'D_', 'E_', 'gi', 'go', 'x8']) and v.x == 1:
            File.write(v.varName)
            File.write('\n')
    File.close()
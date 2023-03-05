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

def Linear_Constr_SPIX(Step):
    obj1 = []
    obj2 = []
    obj3 = []
    Constr = []
    bin_var = []
    general_var = []

    X = [[['X_s' + str(i) + '_' + str(j) + '_' + str(k) for k in range(64)] for j in range(4)] for i in range(Step+1)]
    for i in range(Step+1):
        for j in range(4):
            for k in range(64):
                bin_var.append(X[i][j][k])

    #额外的限制条件
    '''
    SR1 = []
    SR2 = []
    for i in range(32):
        SR1.append(X[0][1][i])
        SR1.append(X[0][3][i])
        SR2.append(X[Step][1][i])
        SR2.append(X[Step][3][i])
    Constr.append(' + '.join(SR1) + ' >= 1')
    Constr.append(' + '.join(SR2) + ' >= 1')

    SC1 = []
    SC2 = []
    for i in range(32, 64):
        SC1.append(X[0][1][i])
        SC1.append(X[0][3][i])
        SC2.append(X[Step][1][i])
        SC2.append(X[Step][3][i])
    for i in range(64):
        SC1.append(X[0][0][i])
        SC1.append(X[0][2][i])
        SC2.append(X[Step][0][i])
        SC2.append(X[Step][2][i])
    Constr.append(' + '.join(SC1) + ' = 0')
    Constr.append(' + '.join(SC2) + ' = 0')
    '''

    #2
    IN = []
    OUT = []
    for i in range(4):
        for j in range(64):
            IN.append(X[0][i][j])
            OUT.append(X[Step][i][j])

    Constr.append(' + '.join(IN) + ' >= 1')
    Constr.append(' + '.join(OUT) + ' >= 1')


    gim = [[[['gim_s' + str(i) + '_' + str(2*j+1) + '_r' + str(k) + '_' + str(l) for l in range(32)] for k in range(8)] for j in range(2)] for i in range(Step)]
    gom = [[[['gom_s' + str(i) + '_' + str(2*j+1) + '_r' + str(k) + '_' + str(l) for l in range(32)] for k in range(8)] for j in range(2)] for i in range(Step)]
    x8 = [[['x8_s' + str(i) + '_' + str(2*j+1) + '_' + str(l) for l in range(32)] for j in range(2)] for i in range(Step)]

    for i in range(Step):
        for j in range(2):
            for k in range(8):
                for l in range(32):
                    bin_var.append(gim[i][j][k][l])
                    bin_var.append(gom[i][j][k][l])

    for i in range(Step):
        for j in range(2):
            for k in range(32):
                bin_var.append(x8[i][j][k])

    #SPIX线性刻画
    for i in range(Step):

        Sout = [X[i+1][3], X[i+1][0], X[i+1][1], X[i+1][2]]

        #对分支的刻画：
        dummy1 = ['dummy_s' + str(i) + '_1' + '_' + str(j) for j in range(64)]
        dummy3 = ['dummy_s' + str(i) + '_3' + '_' + str(j) for j in range(64)]
        for j in range(64):
            general_var.append(dummy1[j])
            general_var.append(dummy3[j])

        branch1 = gom[i][0][7] + x8[i][0]
        branch3 = gom[i][1][7] + x8[i][1]

        for j in range(64):
            Constr = Constr + Linear_Constr_branch([branch1[j], Sout[0][j], Sout[1][j]], dummy1[j])
            Constr = Constr + Linear_Constr_branch([branch3[j], Sout[2][j], Sout[3][j]], dummy3[j])

        #异或刻画
        for j in range(64):
            Constr = Constr + [X[i][0][j] + ' - ' + Sout[0][j] + ' = 0']
            Constr = Constr + [X[i][2][j] + ' - ' + Sout[2][j] + ' = 0']

        #SB-64刻画
        for j in range(2):
            for k in range(8):
                #刻画环型与门组合
                Gim = [None for l in range(32)]
                Gom = [None for l in range(32)]
                for l in range(32):
                    Gim[l] = gim[i][j][k][(l*5)%32]
                    Gom[l] = gom[i][j][k][(l*5)%32]

                L1 = ['L1_s' + str(i) + '_' + str(j) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                L2 = ['L2_s' + str(i) + '_' + str(j) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                Beta = ['Beta_s' + str(i) + '_' + str(j) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                beta = 'beta_s' + str(i) + '_' + str(j) + '_r' + str(k)
                Dummyg = ['Dummyg_s' + str(i) + '_' + str(j) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                dummyg1 = 'dummyg1_s' + str(i) + '_' + str(j) + '_r' + str(k)
                dummyg2 = 'dummyg2_s' + str(i) + '_' + str(j) + '_r' + str(k)

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
                dummysr = ['dummysr_s' + str(i) + '_' + str(j) + '_r' + str(k) + '_' + str(l) for l in range(32)]
                general_var = general_var + dummysr

                branch2 = gim[i][j][k]
                branch3 = [None for l in range(32)]
                for l in range(32):
                    branch3[l] = gom[i][j][k][(l-1)%32]

                if k == 0:
                    branch1 = X[i][2*j+1][0:32]
                else:
                    branch1 = gom[i][j][k-1]
                if k == 7:
                    branch4 = x8[i][j]
                else:
                    branch4 = gom[i][j][k+1]

                for l in range(32):
                    Constr = Constr + Linear_Constr_branch([branch1[l], branch2[l], branch3[l], branch4[l]], dummysr[l])

                #刻画SB-64轮函数中的相等
                for l in range(32):
                    Constr = Constr + [X[i][2*j+1][l+32] + ' - ' + gom[i][j][0][l] + ' = 0']

    return [obj1, obj2, obj3, Constr, bin_var, general_var]

def SPIXLinearmodel_Generate(Step):
    [obj1, obj2, obj3, Constr, bin_var, general_var] = Linear_Constr_SPIX(Step)

    filename = 'SPIXLinearmodel_S' + str(Step) + '.lp'
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

if __name__ == '__main__':
    Step = 8
    SPIXLinearmodel_Generate(Step)
    modelname = 'SPIXLinearmodel_S' + str(Step)
    m = read(modelname)
    m.setParam(GRB.Param.MIPFocus, 2)
    m.optimize()
    vars = m.getVars()
    File = open('Result_SPIXLinear_S' + str(Step) + '.txt', 'w')
    for v in vars:
        File.write(v.varName)
        File.write(':')
        File.write(str(v.x))
        File.write('\n')
    File.close()
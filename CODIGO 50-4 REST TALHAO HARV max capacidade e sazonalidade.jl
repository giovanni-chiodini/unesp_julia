#import Pkg
#Pkg.add("Cbc")
#Pkg.add("JuMP")
#Pkg.add("Statistics")

using Cbc
using JuMP
using Statistics
modelo = Model(Cbc.Optimizer)

T = [16;17;18;19;20;21;22;23;24]#months in which the claims occur
D = [17500;11200;12845;7000;24500;11200;31500;27230;18500] #demandas relacionadas a cada mês
P = [180;190;170;170;160;175;185;170;180;190;
     170;170;160;175;185;170;180;190;170;170;
     160;175;185;170;180;
     240;220;230;250;240;210;220;210;
     240;220;230;250;240;210;220;210;
     240;220;230;250;240;210;220;210;240] #produtividade dos 16 plots
L = [20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;
     20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;
     20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;20.28;] #tamanho de cada plot L
tl = [12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12
      18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18;18] #tempo de colheita da cana
to = [9;9;10;10;9;10;10;9;9;10;10;9;10;10;9;9;10;10;9;10;10;9;9;10;10;
      4;1;2;3;4;2;2;3;1;4;1;2;3;4;2;2;3;1;4;1;2;3;4;2;2] #mês em que ocorreu plantio
J  = [8;13;23;6] #numero de plots em cada fazenda
JF = [(1;2;3;4;5;6;7;8);(9;10;11;12;13;14;15;16;17;18;19;20;21);(22;23;24;25;
        26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44);(45;46;47;48;49;50)]
sazon = [98.78;96.54;98.09;102.11;103.05;102.6;101.89;101.86;100.56]

# as 4 linhas abaixo servem para normalização do problema
v_max =  sum(maximum(T) .- (tl + to)) # calcula o limite superior do objetivo (1)
v_min =  sum(minimum(T) .- (tl + to)) # calcula o limite inferior do objetivo (1)
v_maxt =  maximum((maximum(T) .- (tl + to))) # calcula o limite superior do objetivo (2)
v_mint =  minimum((minimum(T) .- (tl + to))) # calcula o limite inferior do objetivo (2)

m = length(D)
k = length(L)
F = length(J)
G = length(J)
α = 0.5
Cm = 1262.62 #capacidade de colheita da maquina em ton/dia

@variable(modelo, x[i =1:m, j = 1:k], Bin) # variavel binária 1 se há colheita no mês i no plot j
@variable(modelo, y[i =1:m, f = 1:F], Bin) # variavel binária 1 se há colheita no mês i na fazenda f
@variable(modelo, N[i =1:m]) # variavel relacionada colheita nas fazendas em cada mês i
@variable(modelo, t[j = 1:k]) #tempo total colheita
@variable(modelo, d[j = 1:k] >= 0) # desvio de d+ i1/i2
@variable(modelo, dl[j = 1:k] >= 0) # desvio d- de i1/i2
@variable(modelo, θ) #desvio máximo

@constraint(modelo, [j = 1:k], t[j] - to[j] - tl[j] - d[j] + dl[j] == 0); # restrição (5)
@constraint(modelo, [j = 1:k], t[j] == sum(T[i]*x[i,j] for i = 1:m)); #restrição (6)
@constraint(modelo, [j = 1:k], sum(x[i,j] for i = 1:m) == 1); #restrição (7)
@constraint(modelo, [i = 1:m], sum(P[j]*L[j]*x[i,j] for j = 1:k) >= D[i]); #restrição (8)
@constraint(modelo, [i = 1:m, j = 1:k], d[j]+dl[j] <= θ); #restrição (9)
@constraint(modelo, [j = [JF[f] for f = 1:F], f=1:F, i=1:m], x[i,j] <= y[i,f]); #restrição (10)
@constraint(modelo, [i = 1:m], N[i] == sum(y[i,f] for f = 1:F)); #restrição (11)
@constraint(modelo, [i = 1:m], N[i] <= G); #restrição (12)
@constraint(modelo, [i=  1:m], sum(P[j]*L[j]*x[i,j] for j = 1:k) <= Cm*30) # restrição número de talhão colhido.
@constraint(modelo, [i=  1:m], sum( d[j] + dl[j] for j = 1:k) <= 100) # restrição soma dos desvios
@constraint(modelo, [i=  1:m], θ <= 3) # restrição desvio maximo

@objective(modelo, Max, sum(P[j]*L[j]*x[i, j]*(sazon[i]/100) for i in 1:m, j in 1:k))

set_time_limit_sec(modelo, 10.0)

optimize!(modelo)

println("Função objetivo: ", objective_value(modelo))

println("d[j]")
for j=1:k
    println("d[$j] = ", value(d[j]))
end

println("t[j]")
for j=1:k
    println("t[$j] = ", value(t[j]))
end


@time(modelo)
#---------------------------------------------------------------------------------
#Pegando o valor do array de dispersão observadas.
disp = Vector{Float64}()
for j=1:k
append!(disp, value(d[j]))
end

#Pegando o valor do array de número de fazendas observadas por mês
nfazendas = Vector{Float64}()
for i=1:m
append!(nfazendas, value(N[i]))
end

println("
Area total    ", sum(L))
println("Média desvio  ", mean(disp))
println("Desvio maximo ", maximum(disp))
println("% Plots com desvio maior que dois   ", length(disp[disp .> 2.0])/length(disp)*100)
println("Soma desvio   ", sum(disp))
println("Média fazendas observadas por mes   ", mean(nfazendas))
println("Valor da soma dos desvios teta   ", value(θ))
println("Restricao Soma desvios    ", value(sum( d[j] + dl[j] for j = 1:k)))
println("HOY:   ",value(Cm))

#=
for i in 1:m
    produtividade = sum(P[j]*L[j]*x[i,j] for j = 1:k)
    println("Numero de talhoes colhidos no mes $i: $produtividade")
end
=#

#=
for i in 1:m
    num_plots_harvested = sum(value(x[i, j]) for j in 1:k)
    println("Numero de talhoes colhidos no mes $i: $num_plots_harvested")
end
=#


for i in 1:m
    println("
Talhoes colhidos no mes $i:")
    for j in 1:k
        if value(x[i, j]) > 0.5  # Verifica se o plot j foi colhido no mes i
            println("- $j")
        end
    end
end

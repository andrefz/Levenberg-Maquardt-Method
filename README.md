# Método de Levenberg-Marquardt

Este repositório contém implementações em Julia 1.8.x do Método Levenberg 
Marquardt para problema de quadrados mínimos.

[Relatório com testes de desempenho][relatorio-lm]

[relatorio-lm]: https://github.com/andrefz/Levenberg-Marquardt-Method/blob/main/tex/Levenberg_Maquardt_Method.pdf


## Reprodutibilidade

Para reproduzir os experimentos relatados, basta executar o arquivo
`Levenberg_Marquardt.ipynb` em um ambiente `IJulia`.

O arquivo `algorithm.jl` disponibiliza uma implementação com diferenciação
automática e outra cuja derivada do modelo é parâmetro.   

## Autores
Os autores do repósitorio, listados em ordem alfabética, são:

* [@andrefz](https://github.com/andrefz) - _André F. Zanella_
* [@juliaguizardi](https://github.com/JuliaGuizardi) - _Julia D. Guizardi_

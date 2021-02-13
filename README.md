# Filogenia-para-iniciantes
 Este tutorial foi criado por Olivier De Clerk (Ugent). Atualizado e Organizado por Daniel Moura


### Árvores: estruturas básicas
R oferece uma ampla gama de pacotes para trabalhar com filogenias e dados moleculares. Uma visão geral quase abrangente pode ser encontrada em - https://cran.r-project.org/web/views/Phylogenetics.html. Para ser mais didático, vamos nos concentrar em algumas funções básicas incluídas nos pacotes mais comumente usados.

Começaremos com “ape”, que abreviação do inglês Análises da Filogenética e Evolução (*Analyses of Phylogenetics and Evolution*). “ape” faz muitas coisas diferentes, principalmente criar árvores filogenéticas. Informações sobre “ape” podem ser encontradas no livro: Paradis, E. 2011. Analysis of Phylogenetics and Evolution with R. Springer Science & Business Media. Os livros USE-R são publicados pela Springer(http://www.springer.com/series/6991).

Para começar, carregue a biblioteca e simule e plote uma filogenia, que a "ape" pode fazer em diferentes modelos (Lembre de instalar ela antes).
```
# para instalar digite 
install.packages("ape")

# só precisa ser feita uma única vez na máquina

# para carregar a biblioteca "ape"
library(ape)
```
e provavelmente receberá essa resposta:

> Warning: package 'ape' was built under R version 3.1.3

Em seguida execute o seguite código:
```
tree <- rtree(n = 20)
plot(tree, edge.width = 2)
```
a função rtree irá gerar uma árvore aleátoria, onde n=20 indica o numero de taxons desejado. Na função plot irá reproduzir a árvore visualmente e com a opção edge.width é possível selecionar a largura das linhas.

![image](https://user-images.githubusercontent.com/23074735/89691809-4a5b6900-d8e0-11ea-8218-423dafa12dc8.png)

a função rtree irá gerar uma árvore aleátoria, onde n=20 indica o número de taxons desejado. Na função "plot" irá reproduzir a árvore visualmente e com a opção "edge.width" é possível selecionar a largura das linhas. Tente utilizar tamanhos de taxons e de larguras diferentes para compreender o código.

Isso pode ser uma Simulação ou estimação de uma filogenia, ou leitura um de um arquivo de entrada, que é uma lista da classe "phylo". Uma lista é apenas um tipo de objeto personalizável que pode combinar diferentes objetos de diferentes tipos. Por exemplo, uma lista pode ter um vetor de números reais (com a classe "numérica") como seu primeiro elemento; e então um vetor de strings (com classe "caractere") como seu segundo elemento; e assim por diante.

Um objeto da classe "phylo" possui pelo menos 4 partes. Eles normalmente estão ocultos, por exemplo, apenas digitar o nome do seu objeto "phylo" não fornece a estrutura na memória, como acontece com muitos objetos R. Para chegar ao entendimento de como um objeto da classe "phylo" codifica filogenias, vamos usar um caso simples: uma árvore com 5 pontas e sem comprimentos de borda.

```
tree <- read.tree(text = "(((A,B),(C,D)),E);")
plot(tree, edge.width = 2, type = "cladogram",  font = 4)
```

Dento de "text" se coloca os parênteses para indicar um agrupamento dos grupos, ou seja, é uma árvore em forma de script.

A árvore é composta de arestas (ramificação, no inglês, branches), que conectam os nós (ponto de ramificaçã) e a ponta. Você pode inspecionar as estruturas e visualizá-las usando os comandos a seguir.

```
tree$edge
```

A aresta da matriz contém o número do nó inicial e final para todos os nós e pontas da árvore. Por convenção, as pontas da árvore são numeradas de 1 a n para n pontas; e os nós são numerados de n + 1 a n + m para m nós. m = n - 1 para uma árvore totalmente bifurcada. Isso é apenas para manter o controle de quais nós são internos e quais são folhas.

```
tree$tip.label
```

O vetor tip.label contém os rótulos de todas as pontas na árvore. A ordem de tip.label é a ordem das pontas numeradas de 1 a n na borda.

```
plot(tree, type = "phylogram", edge.width = 2)
tiplabels()
nodelabels()
```
Dê uma olhada na numeração das pontas(tips) e dos nós(node).

Podemos ver que todas as informações de relacão entre os táxons na árvore estão contidas nos nós inicial e final de cada aresta. As arestas que compartilham um número de nó inicial comum são descendentes de um ancestral comum imediato, etc.

A árvore também pode ter outros componentes, os mais comuns dos quais são edge.length (um vetor da classe "numérico" contendo todos os comprimentos das arestas da árvore na mesma ordem das linhas na aresta; e root.edge, um valor numérico que fornece o comprimento da borda da raiz, se houver. Além disso, outros elementos e atributos podem ser adicionados para tipos especiais de árvores filogenéticas.
                                                                                                           
É importante ter em mente que esta não é a única maneira de armazenar uma filogenia na memória do computador, nem mesmo a única maneira de armazenar uma árvore em R! No entanto, a vantagem de armazenar árvores da mesma maneira em pacotes diferentes é que isso melhora a interoperabilidade de pacotes - uma maneira elegante de dizer que os desenvolvedores podem se concentrar no desenvolvimento de métodos que se baseiam na funcionalidade existente, em vez de reinventar a roda novamente para cada projeto.
                                                                                                           
Por exemplo, uma das muitas razões para construir sobre as estruturas desenvolvidas para o pacote "ape", é que existem muitas funções utilitárias diferentes em "ape" e os outros pacotes R para filogenética (por exemplo, "phytools"), para leitura, escrita e manipulação de árvores filogenéticas armazenadas na memória usando essa estrutura. Há um grande número de funções de "utilidade" desse tipo para filogenias, então vou apenas revisar algumas das mais importantes e úteis.


### Escrevendo e lendo árvores filogenéticas

Para ilustrar os comandos básicos relacionados à leitura e escrita de árvores, usaremos o pacote "phytools". Lembre de instalar o pacote antes de executar.

```
# para instalar
install.packages("phytools")

# só precisa ser feita uma única vez na máquina

# para carregar a biblioteca "phytools"
library(phytools)
```
Em seguida execute o comando para carregar uma árvore aleatória.

```
tree2 <- rtree(n = 18)
plot(tree2, edge.width = 2)
```
Agora é preciso criar um arquivo que salve essa informação da árvore gerada.

```
writeNexus(tree2,"minha_primeira_arvore.nex")
```
Você pode visualizar o arquivo nexo usando o comando cat, que lê e imprime as linhas do arquivo texto.
```
cat(readLines("minha_primeira_arvore.nex"), sep = "\n")
```
Também podemos extrair clados e excluir pontas.

Para simular uma árvore de birth-death usando "phytools". b e d são as taxas de especiação e extinção (b = birth e d= death),e n é o número de taxon.
```
tree <- pbtree(b = 1, d =0.2 , n = 40)
plotTree(tree)
nodelabels()
```

Agora podemos extrair o clado descendente de um nó específico (por exemplo, o nó 66, ou qualquer outro nó que produz descendentes suficientes e vamos salvá-lo como uma árvore separada).

```
clado80 = extract.clade(tree,80)
plotTree(clado80)
```

Então podemos salvar somente esse modelo do nó 66.

```
writeNexus(clado80, "node80.nex")
```

Vamos tirar 10 pontas da árvore original (aleatoriamente).

```
dtips <- sample(tree$tip.label, 10)
dt <- drop.tip(tree, dtips)
plotTree(dt)
```
A função drop.tip irá remover conforme indicado na segunda posição (após a vírgula). Nesse caso, representado pela variavável "dt" que foi criada com pontas aleatórias.

Também podemos descartar todas as pontas que foram extintas antes do presente. Este método abaixo é legal, porém não é o único.

```
et <- fancyTree(tree, type = "droptip", tip = getExtinct(tree), cex = 0.7)
```

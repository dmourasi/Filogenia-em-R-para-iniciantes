# Filogenia-para-iniciantes
 Este tutorial foi criado por Olivier De Clerk (Ugent). Atualizado e Organizado por Daniel Moura


### Árvores: estruturas básicas
R oferece uma ampla gama de pacotes para trabalhar com filogenias e dados moleculares. Uma visão geral quase abrangente pode ser encontrada em - https://cran.r-project.org/web/views/Phylogenetics.html. Para ser mais didático, vamos nos concentrar em algumas funções básicas incluídas nos pacotes mais comumente usados.

Começaremos com “ape”, que abreviação do inglês Análises da Filogenética e Evolução (*Analyses of Phylogenetics and Evolution*). “ape” faz muitas coisas diferentes, principalmente criar árvores filogenéticas. Informações sobre “ape” podem ser encontradas no livro: Paradis, E. 2011. Analysis of Phylogenetics and Evolution with R. Springer Science & Business Media. Os livros USE-R são publicados pela Springer(http://www.springer.com/series/6991).

Para começar, carregue a biblioteca e simule e plote uma filogenia, que a "ape" pode fazer em diferentes modelos (Lembre de instalar ela antes).
```
#carregar a biblioteca "ape"
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

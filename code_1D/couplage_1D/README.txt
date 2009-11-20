Dans ce repertoire sont places les binaires ALOHA, 
compilés à partir du Fortran 77 ou du Fortan 90.

Le programme matlab fait la distinction entre les binaires : 
pour cela, il regarde sur quel plateforme (32bits, 64bits, etc)
il se trouve, et appelle le binaire correspondant.

Pour se faire, les noms de fichiers doivent correspondre aux regles
suivantes : 

coupl_plasma_versionX_ARCHITECTURE

ou le X correspond à la version du code (3,6...)
et ARCHITECTURE correspond à l'architecture sur laquelle
le code a été compilé : 
- glnxa64 (nashira) 
- alpha (deneb)
- glnx86 (x86)

Bootstrap: docker
From: thomaschln/r-devtools

%post
  apt install -y git
  cd /opt
  git clone https://github.com/hmgu-itg/man_qq_annotate.git -b package
  cd man_qq_annotate
  R -e 'library(devtools); install()'

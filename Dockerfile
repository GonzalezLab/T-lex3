FROM debian:8-slim

# Install & update required system dependencies
RUN apt-get update && \
  apt-get dist-upgrade -y -qqq && \
  apt-get autoremove -y -qqq && \
  apt-get install -y -qqq perl bash wget unzip bzip2 libkrb5-3 libpng12-0 && \
  apt-get clean -y -qqq && \
  rm -rf /var/lib/apt/lists/*

#
# INSTALL RMBLAST & BLAST+
#
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz && \
  wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
  tar -xzf ncbi-rmblastn* && \
  tar -xzf ncbi-blast* && \
  mv ncbi*/bin/* bin && \
  rm -rf ncbi*

#
# INSTALL TRF (FOR REPEATMASKER)
#
RUN cd /usr/local/bin && \
  wget http://tandem.bu.edu/trf/downloads/trf407b.linux64 && \
  mv trf*.linux64 trf && \
  chmod +x trf

#
# INSTALL REPEATMASKER
#
RUN cd /usr/local && \
  wget http://www.repeatmasker.org/RepeatMasker-open-4-0-9-p2.tar.gz && \
  tar -xzf RepeatMasker-open-4-0-9-p2.tar.gz && \
  rm -rf RepeatMasker-open-4-0-9-p2.tar.gz && \
  perl -0p -e 's/\/usr\/local\/hmmer/\/usr\/bin/g;' \
  -e 's/\/usr\/local\/rmblast/\/bin/g;' \
  -e 's/DEFAULT_SEARCH_ENGINE = "crossmatch"/DEFAULT_SEARCH_ENGINE = "ncbi"/g;' \
  -e 's/TRF_PRGM = ""/TRF_PRGM = "\/usr\/local\/bin\/trf"/g;' RepeatMasker/RepeatMaskerConfig.tmpl > RepeatMasker/RepeatMaskerConfig.pm

# Fix RepeatMasker's strange shebang lines
RUN cd /usr/local/RepeatMasker \
  && perl -i -0pe 's/^#\!.*perl.*/#\!\/usr\/bin\/env perl/g' \
  RepeatMasker \
  DateRepeats \
  ProcessRepeats \
  RepeatProteinMask \
  DupMasker \
  util/queryRepeatDatabase.pl \
  util/queryTaxonomyDatabase.pl \
  util/rmOutToGFF3.pl \
  util/rmToUCSCTables.pl

# Allow us to find RepeatMasker binaries
ENV PATH=$PATH:/usr/local/RepeatMasker

#
# REPBASE LIBRARIES
#
COPY ./RepBase.tar.gz /usr/local/RepeatMasker/
RUN cd /usr/local/RepeatMasker && \
  tar -xzf RepBase.tar.gz && \
  rm -rf RepBase.tar.gz


#
# INSTALL SHRIMP 2
#
RUN cd /usr/local && \
  wget http://compbio.cs.toronto.edu/shrimp/releases/SHRiMP_2_2_3.lx26.x86_64.tar.gz && \
  tar -xzf SHRiMP_2_2_3.lx26.x86_64.tar.gz && \
  rm -rf SHRiMP_2_2_3.lx26.x86_64.tar.gz

ENV SHRIMP_PATH /usr/local/SHRiMP_2_2_3
ENV PATH $PATH:$SHRIMP_PATH/bin:$SHRIMP_PATH/utils

#
# INSTALL BLAT
#
RUN cd /usr/local && \
  mkdir blat36 && \
  wget https://hgwdev.gi.ucsc.edu/~kent/exe/linux/blatSuite.36.zip && \
  unzip blatSuite.36.zip -d /usr/local/blat36 && \
  rm -rf blatSuite.36.zip

ENV BLAT_PATH /usr/local/blat36
ENV PATH $PATH:$BLAT_PATH

#
# INSTALL SAMTOOLS & BCFTOOLS & BWA
#
RUN apt-get update -q && \
  apt-get install -y -qqq build-essential libncurses-dev libbz2-dev liblzma-dev zlib1g-dev && \
  cd /tmp && \
  wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
  tar -xjf samtools-1.8.tar.bz2 && \
  cd samtools-1.8 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd /tmp && \
  wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2 && \
  tar -xjf bcftools-1.8.tar.bz2 && \
  cd bcftools-1.8 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  rm -rf /var/lib/apt/lists/* && \
  cd /tmp && \
  wget https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2 && \
  tar -xjf bwa-0.7.17.tar.bz2 && \
  cd bwa-0.7.17 && \
  make && \
  cp bwa /usr/local/bin/ && \
  apt-get remove -y -qqq build-essential libncurses-dev libbz2-dev liblzma-dev && \
  apt-get autoremove -y -qqq && \
  apt-get clean -y -qqq && \
  cd /tmp && \
  rm -rf /tmp/*

#
# INSTALL FASTX TOOLKIT & DEPENDENCIES
#
RUN apt-get update -q && \
  apt-get install -y -qqq build-essential pkg-config autoconf && \
  cd /tmp && \
  wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz && \
  tar -xzf libgtextutils-0.7.tar.gz && \
  cd libgtextutils-0.7 && \
  ./reconf && \
  ./configure && \
  make && \
  make install && \
  cd /tmp && \
  wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2 && \
  tar -xjf fastx_toolkit-0.0.14.tar.bz2 && \
  cd fastx_toolkit-0.0.14 && \
  ./reconf && \
  ./configure && \
  make && \
  make install && \
  apt-get remove -y -qqq build-essential pkg-config autoconf && \
  apt-get autoremove -y -qqq && \
  apt-get clean -y -qqq && \
  cd /tmp && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /tmp/*

#EMBOSSFASTX TOOLKIT & DEPENDENCIES
#
RUN apt-get update -q && \
  apt-get install -y -qqq emboss && \
  apt-get clean -y -qqq && \
  cd /tmp && \
  rm -rf /var/lib/apt/lists/*

#
# INSTALL FASTAGREP
#
RUN cd /tmp && \
  wget http://bioinfo.ut.ee/download/dl.php?file=28 -O fastagrep.tar.gz && \
  tar -xzf fastagrep.tar.gz && \
  mv ./fastagrep_v2.0_64bit_linux_2_6 /usr/local/bin/fastagrep && \
  rm -rf fastagrep.tar.gz

#
# INSTALL TLEX
#

RUN mkdir -p /data
VOLUME [ "/data" ]
WORKDIR /data

COPY ./tlex-open-v3.0.pl /usr/bin/tlex
RUN chmod +x /usr/bin/tlex

ENTRYPOINT [ "/usr/bin/tlex" ]

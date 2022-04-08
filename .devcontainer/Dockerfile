# [Choice] Ubuntu version: bionic, focal
FROM mcr.microsoft.com/vscode/devcontainers/base:ubuntu-20.04

# [Optional] Uncomment this section to install additional OS packages.
#RUN apt-get update && export DEBIAN_FRONTEND=noninteractive
#RUN apt-get -y install git
# RUN bash -i -c 'pip3 install --user --no-use-pep517 -e .[all]'
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install make && \
    apt-get install --no-install-recommends -y python3.8 python3-pip python3.8-dev pylint

# [Optional] Uncomment this line to install global node packages.
# RUN su vscode -c "source /usr/local/share/nvm/nvm.sh && npm install -g <your-package-here>" 2>&1
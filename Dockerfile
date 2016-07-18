FROM centos:centos7

MAINTAINER Peyton McNully HudsonAlpha Institute  <pmcnully@hudsonalpha.org>

RUN yum -y update; yum clean all
RUN yum -y install epel-release; yum clean all
RUN yum -y install python-pip; yum clean all
RUN yum -y install R; yum clean all
RUN yum -y install curl; yum clean all

ENV org=HudsonAlpha repo=aRNAPipe api=api.github.com

RUN LATEST_RELEASE=$(curl -s https://${api}/repos/${org}/{repo}/releases \
 | grep browser_download_url \
 | head -n 1 | cut -d '"' -f 4)

ADD ${LATEST_RELEASE} /src

RUN cd /src; pip install -r requirements.txt

EXPOSE 8080

CMD ["python", "/src/index.py"]

FROM postgres:15 

ARG PG_POSTGRES_PWD=smadi
ARG DBUSER=user1
ARG DBUSER_PWD=smadi
ARG DBNAME=sampledb
ARG DB_DUMP_FILE=input/init.sql

ENV POSTGRES_DB launchpad
ENV POSTGRES_USER postgres
ENV POSTGRES_PASSWORD ${PG_POSTGRES_PWD}
ENV PGDATA /pgdata

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN apt-get update
RUN apt-get dist-upgrade -y
RUN apt-get autoremove -y
RUN apt-get install tzdata locales -y
RUN ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen
RUN dpkg-reconfigure tzdata
RUN locale-gen



RUN apt-get update && apt-get install -y python3-pip && pip3 install --break-system-packages sympy && pip3 install --break-system-packages numpy && pip3 install --break-system-packages scipy 
RUN   apt-get install -y postgresql-plpython3-15 

RUN  apt-get install -y dos2unix
COPY jdbc/wait-for-pg-is-ready.sh /tmp/wait-for-pg-is-ready.sh
RUN dos2unix /tmp/wait-for-pg-is-ready.sh && apt-get --purge remove -y dos2unix && rm -rf /var/lib/apt/lists/*


COPY ${DB_DUMP_FILE} /tmp/init.sql

RUN set -e && \
    nohup bash -c "docker-entrypoint.sh postgres &"  && \
	  /tmp/wait-for-pg-is-ready.sh && \
    psql -U postgres -c "CREATE USER ${DBUSER} WITH SUPERUSER CREATEDB CREATEROLE ENCRYPTED PASSWORD '${DBUSER_PWD}';" && \
    psql -U ${DBUSER} -d ${POSTGRES_DB} -c "CREATE DATABASE ${DBNAME} TEMPLATE template0;" && \
    psql -U ${DBUSER} -d ${DBNAME} < /tmp/init.sql && \
    psql -U ${DBUSER} -d ${DBNAME} -c "REFRESH MATERIALIZED VIEW public.concr_parameter;"

HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
  CMD pg_isready -U postgres -d launchpad
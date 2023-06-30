FROM python:3.10

LABEL name=a290dash \
    maintainer="Jooa <j.hooli@dkfz.de>" \
    description="Interactive single-cell visualisation" \
    version="0.1.0"

WORKDIR /usr/src/a290dash

ADD app.py .
ADD assets assets/
ADD requirements.txt .

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV PATH "$PATH:/usr/local/bin"
ENV DASH_DATA_DIR "/data"
VOLUME /data

RUN pip install -r requirements.txt


CMD [ "/usr/local/bin/gunicorn", "app:SERVER",  "-b", ":8000" ]

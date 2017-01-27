<template lang='pug'>
div
  h1 {{title}}
  slot(name='head')
  form.center-block
    honeypot(v-model='honeypot')
    slot(name='form')
    .btn.btn-lg.btn-default.center-block.submit-button(v-if='submitButton && (!polling.inProgress)' @click='submit()') {{submitButtonText }}

  .polling(v-if='polling.inProgress')
    spinner(:color="spinner.color", :size="spinner.size")
    .polling-message {{polling.message}}
    el-steps(:space='100', :active='progressStage')
      el-step(icon='upload')
      el-step(icon='setting')
      el-step.last-step(icon='message')

  el-alert(v-if='result.error', :title="result.error", type="warning")
  el-alert(v-if='requestError', :title="requestError", type="error")


  .results(v-if='!polling.inProgress')
    .results-summary(v-if='result.preview', v-html="result.preview.html")
    downloadbutton(:filedata='result.file' v-if='result.file')
</template>

<script>
import honeypot from '../../components/widgets/Honeypot'
import downloadbutton from '../../components/widgets/DownloadButton'
import spinner from 'vue-spinner/src/PulseLoader'
import tools from '../../tools.js'

export default {
  mixins: [tools.valueTracker],
  props: {
    title: {default: 'Title here'},
    submitButton: {default: true},
    submitButtonText: {default: 'Submit'},
    queryUrl: {default: '/'},
    value: {default: () => ({})},
    backendUrl: {default: ''}
  },
  data: function () {
    return {
      fields: {},
      honeypot: '',
      polling: {
        inProgress: false,
        status: 'queued',
        message: '',
        data: ''
      },
      result: {},
      requestError: '',
      spinner: {
        size: '12px',
        color: '#6da5ff'
      }
    }
  },
  components: {
    honeypot,
    downloadbutton,
    spinner
  },
  watch: {
    fields: {
      handler: function (val) {
        console.log(val)
        this.$emit('input', val)
      },
      deep: true
    },
    value: {
      handler: function (val) {
        this.fields = val
      },
      deep: true
    }
  },
  computed: {
    fulldata: function () {
      return {
        data: this.fields,
        honeypot: this.honeypot
      }
    },
    progressStage: function () {
      return ['queued', 'started', 'finished'].indexOf(this.polling.status) + 1
    }
  },
  methods: {
    submit: function () {
      console.log('fields', this.fields)
      console.log('sc2', this.$http, this.backendUrl)
      console.log('sc3', 'http://localhost:8000/' + this.backendUrl)
      this.$http.post(
        'http://localhost:8000/' + this.backendUrl,
        this.fields
      ).then(function (response) {
        // SUCCESS
        console.log('res', response)
        console.log('jobstart response', response)
        if (response.success === false) {
          return
        }
        this.polling.inProgress = true
        var jobPoller = setInterval(function () {
          this.$http.post(
            'http://localhost:8000/poll',
            {job_id: response.body.job_id}
          ).then(function (pollResponse) {
            console.log('poll response', pollResponse)
            // clearInterval(jobPoller)
            this.request_error = pollResponse.error
            var data = pollResponse.body
            this.polling.status = data.status
            if (data.success === false) {
              clearInterval(jobPoller)
              this.polling.inProgress = false
            } else if (data.status === 'queued') {
              this.polling.message = 'Job pending...'
            } else if (data.status === 'started') {
              this.polling.message = data.progress.message
              this.polling.data = data.progress.data
            } else if (data.status === 'finished') {
              this.polling.message = 'Finished ! sending now'
              this.result = data.result
              this.polling.inProgress = false
              clearInterval(jobPoller)
            }
            this.polling.message = data.progress.message
            console.log('polmess', this.polling.message)
            // clearInterval(jobPoller)
          }, function (response) {
            // FAILURE OF THE POLLING
            this.polling.inProgress = false
            console.log('polling error res', response)
            this.requestError = 'Sorry, we met an error while polling job status'
            clearInterval(jobPoller)
          })
        }.bind(this), 100)
      }, function (response) {
        // FAILURE OF THE JOB STARTING
        console.log('job starting error res', response)
        if (response.status === 400) {
          this.requestError = 'Bad Request: ' + window.JSON.stringify(response.body)
        }
      })
    }
  }
}
</script>

<style scoped>

.center {
  text-align: center;
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;

}

form {
  width:500px;
  max-width:90%;
  margin-left: auto;
  margin-top: 2em;
}

form h4 {
  text-transform: uppercase;
  color: #000;
  text-align: center;
  position: relative;
  font-size: 18px;
  padding-top: 15px;
  padding-bottom: 15px;
}

.submit-button {
  margin-top: 20px;
  margin-bottom: 60px;
  width: 150px;
}

.el-alert {
  margin: 15px auto 15px;
  width: auto;
  max-width: 500px;
  padding-right: 40px
}

.polling {
  text-align:center;
  margin-top: 50px
}

.last-step {
  width: 40px !important;
}

</style>

<template lang='pug'>
div
  honeypot(v-model='honeypot')
  .btn.btn-lg.btn-default.center-block.submit-button(v-if='submitButton && (!polling.inProgress)' @click='submit()') {{submitButtonText }}
  .polling(v-if='polling.inProgress')
    spinner(:color="spinner.color", :size="spinner.size")
    .polling-message {{polling.message}}
    el-steps(:space='100', :active='progressStage')
      el-step(icon='upload')
      el-step(icon='setting')
      el-step.last-step(icon='message')

  //- el-alert(v-if='result.error', :title="result.error", type="warning")
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
    submitButton: {default: true},
    submitButtonText: {default: 'Submit'},
    value: {default: () => ({})},
    backendUrl: {default: ''},
    backendRoot: {default: 'http://localhost:8000/'},
    form: {default: () => ({})}
  },
  data: function () {
    return {

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
    form: {
      handler: function (val) {
        console.log(val)
        this.$emit('input', val)
      },
      deep: true
    },
    value: {
      handler: function (val) {
        this.form = val
      },
      deep: true
    }
  },
  computed: {
    fulldata: function () {
      return {
        data: this.form,
        honeypot: this.honeypot
      }
    },
    progressStage: function () {
      return ['queued', 'started', 'finished'].indexOf(this.polling.status) + 1
    }
  },
  methods: {
    submit: function () {
      console.log('form', this.form)
      console.log('sc2', this.$http, this.backendUrl)
      this.$http.post(
        this.backendRoot + this.backendUrl,
        this.form
      ).then(function (response) {
        // SUCCESS
        console.log('res', response)
        console.log('jobstart response', response)
        this.startPolling(response.body.job_id)
      }, function (response) {
        // FAILURE OF THE JOB STARTING
        console.log('job starting error res', response)
        this.requestError = 'Failed with status ' + response.status
        if (response.status === 400) {
          this.requestError = 'Bad Request: ' + window.JSON.stringify(response.body)
        }
      })
    },
    startPolling: function (jobId) {
      this.polling.inProgress = true
      var jobPoller = setInterval(function () {
        this.$http.post(
          this.backendRoot + 'poll',
          {job_id: jobId}
        ).then(function (pollResponse) {
          console.log('poll response', pollResponse)
          // clearInterval(jobPoller)
          this.requestError = pollResponse.error
          var data = pollResponse.body
          this.polling.status = data.status
          if (data.success === false) {
            this.requestError = data.error
            this.polling.inProgress = false
            clearInterval(jobPoller)
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
        }, function (pollResponse) {
          // FAILURE OF THE POLLING
          this.polling.inProgress = false
          console.log('polling error res', pollResponse)
          this.requestError = 'Failed with status ' + pollResponse.status
          clearInterval(jobPoller)
        })
      }.bind(this), 100)
    }
  }
}
</script>

<style scoped>

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

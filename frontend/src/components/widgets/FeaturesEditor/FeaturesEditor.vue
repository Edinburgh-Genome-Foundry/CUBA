<template lang='pug'>
.features-editor
  p.label Features minimap
  features-minimap(:featuresData='featuresData', :window='window',
                   :sequenceLength='sequence.length', @drag='updateWindowCenter',
                   @select='selectFeature', @mouseover="hoverFeature")
  p.hoveredFeature
    span(v-if='hoveredFeature') {{features[hoveredFeature].label}}
  p.label Sequence
  sequence-viewer(:sequence='sequence', :window='window')
  features-viewer(:featuresData='featuresData', :window='window',
                  :sequenceLength='sequence.length',
                  @drag='shiftWindowCenter'
                  @select='selectFeature',
                  @mouseover="hoverFeature")

  .editor
    p.label Feature editor
    el-row.location
      el-col(:xs='12', :sm='8', :lg='8')
        p Start
        el-input-number(v-model='editorFeature.start', size='small')
      el-col(:xs='12', :sm='8', :lg='8')
        p End
        el-input-number(v-model='editorFeature.end', size='small')
      el-col(:xs='12', :sm='8', :lg='8')
        p Strand
        el-select(v-model='editorFeature.strand', placeholder='Strand', size='small')
          el-option(label='None', :value='0')
          el-option(label='+', :value='1')
          el-option(label='-', :value='-1')
    p Feature label
    el-input(v-model='editorFeature.label', placeholder='Annotation Label')
    .buttons
      el-button(size='small', icon="plus" @click='newFeature') New
      el-button(v-if='selectedFeature' size='small', icon='edit' @click='updateFeature(selectedFeature.id)') Update
      el-button(v-if='selectedFeature' size='small', icon='delete' @click='deleteFeature(selectedFeature.id)') Delete
</template>
<script>
import sequenceviewer from './SequenceViewer'
import featuresviewer from './FeaturesViewer'
import featuresminimap from './FeaturesMinimap'

export default {
  name: 'features-editor',
  props: {
    sequence: {default: ''},
    value: {default: () => ({})},
    featureColor: {default: () => () => '#20a0ff'}
  },
  data () {
    return {
      lineHeight: 3,
      maxLines: 8,
      // sequence: 'ATGC'.repeat(200),
      clientWidth: 500,
      zoomFactor: 15,
      windowCenter: 0,
      editorFeature: {
        id: null,
        start: 0,
        end: 1,
        strand: 1,
        label: ''
      },
      features: Object.assign({}, this.value),
      hoveredFeature: null
      // featuresMirror: this.features
    }
  },
  components: {
    'sequence-viewer': sequenceviewer,
    'features-viewer': featuresviewer,
    'features-minimap': featuresminimap
  },
  computed: {
    featuresData () {
      console.log('recomputing', this.features)
      var self = this
      return Object.keys(this.features).map(function (featureID) {
        console.log(featureID, self.editorFeature.id)
        var dict = Object.assign({
          level: self.featureLevels[featureID]
        }, self.features[featureID])
        dict['selected'] = (featureID === self.editorFeature.id)
        dict['js_color'] = self.featureColor(dict)
        return dict
      })
    },
    selectedFeature () {
      return this.features[this.editorFeature.id]
    },
    featureLevels () {
      var self = this
      var featureOverlaps = {}
      var featureIDs = Object.keys(this.features)
      featureIDs.map(function (id) {
        featureOverlaps[id] = []
      })
      for (var i = 0; i < featureIDs.length; i++) {
        var id1 = featureIDs[i]
        var f1 = this.features[id1]

        for (var j = i + 1; j < featureIDs.length; j++) {
          var id2 = featureIDs[j]
          var f2 = this.features[id2]
          console.log(f1.start, f1.end, f2.start, f2.end)
          if (((f1.start <= f2.start) && (f2.start < f1.end)) ||
              ((f2.start <= f1.start) && (f1.start < f2.end))) {
            featureOverlaps[id1].push(id2)
            featureOverlaps[id2].push(id1)
          }
        }
      }
      console.log('overlaps', featureOverlaps)
      var levels = {}
      var sortedFeatureIDs = featureIDs.sort(function (id1, id2) {
        var f1 = self.features[id1]
        var f2 = self.features[id2]
        return (f2.end - f2.start) - (f1.end - f1.start)
      })
      sortedFeatureIDs.map(function (id) {
        console.log(id)
        var level = 0
        while (featureOverlaps[id].some((otherID) => (levels[otherID] === level))) {
          level++
        }
        levels[id] = level
      })
      return levels
    },
    window () {
      var windowSize = Math.floor(this.clientWidth / this.zoomFactor)
      windowSize = windowSize + (windowSize % 2)
      var window = {
        start: this.windowCenter - windowSize / 2,
        end: this.windowCenter + windowSize / 2
      }
      if (window.start < 0) {
        window.start = 0
        window.end = windowSize
      } else if (window.end > this.sequence.length) {
        window.end = this.sequence.length
        window.start = this.sequence.length - windowSize
      }
      return window
    }
  },
  methods: {
    handleResize () {
      this.clientWidth = this.$el.parentElement.clientWidth
    },
    updateWindowCenter (center) {
      this.windowCenter = center
    },
    shiftWindowCenter (shift) {
      this.windowCenter += shift
    },
    selectFeature (featureID) {
      console.log('wwwwwww', featureID)
      this.editorFeature = Object.assign({}, this.features[featureID])
    },
    newFeature () {
      var newID = 0
      while (this.features[newID]) {
        newID++
      }
      newID = newID.toString()
      this.$set(this.features, newID, {
        id: newID,
        start: this.editorFeature.start,
        end: this.editorFeature.end,
        strand: this.editorFeature.strand,
        label: this.editorFeature.label
      })
      this.selectFeature(newID)
    },
    updateFeature (featureID) {
      console.log('update', featureID)
      this.features[featureID] = Object.assign({}, this.editorFeature)
    },
    deleteFeature (featureID) {
      this.$delete(this.features, featureID)
    },
    hoverFeature (featureID) {
      console.log(featureID)
      this.hoveredFeature = featureID
    }
  },
  mounted () {
    window.addEventListener('resize', this.handleResize)
    this.handleResize()
  },
  watch: {
    features: {
      deep: true,
      handler (value) {
        this.$emit('change', value)
      }
    },
    value: {
      deep: true,
      handler (value) {
        this.features = value
      }
    }
  }
}
</script>
<!-- Add "scoped" attribute to limit CSS to this component only -->
<style lang='scss' scoped>
.features-editor {
  width: 100%;

  margin-top: 2em;
  border-radius: 5px;
  border: 2px solid #f4faff;
  /*padding: 8px;*/

  p.label {
    font-size: 14px;
    margin-bottom: 0;
    font-weight: bold;
  }
  p.hoveredFeature {
    font-size: 14px;
    font-family: 'Inconsolata';
    text-align: center;
    height: 1em;
    width: 90%;
    margin-left: 5%;
  }
  .editor {
    padding: 5px;
    background-color: #fbfdff;
    .buttons {
      font-size: 16px;
      padding: 3px;
      text-align: center;
    }

    p {
      margin-top: 0;
      font-size: 13px;
      margin-bottom: 0;
      margin-top: 1em;
    }
  }
}
#sequence-vuer {
    height: 300px;
}

</style>
